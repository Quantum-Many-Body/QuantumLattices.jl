module FockPackage

using StaticArrays: SVector
using LinearAlgebra: dot
using Printf: @printf, @sprintf
using ..Spatials: PID, AbstractBond, Point, Bond, decompose
using ..DegreesOfFreedom: IID, Internal, Index, LaTeX, OID, AbstractCompositeOID, latexformat, OIDToTuple, Operator, Operators, Config, Table, oidtype
using ..Terms: wildcard, Subscripts, @subscripts_str, SubID, Coupling, Couplings, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, delta, decimaltostr
using ...Prerequisites.Traits: rawtype
using ...Mathematics.AlgebraOverFields: SimpleID, ID, Element
using ...Mathematics.VectorSpaces: CartesianVectorSpace

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: subscriptsexpr, nonconstrain, otype, couplingcenters, abbr, termfactor
import ...Interfaces: rank, ⊗, ⋅, expand, expand!, permute
import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent

export ANNIHILATION, CREATION, MAJORANA, fdefaultlatex, bdefaultlatex, usualfockindextotuple, nambufockindextotuple
export FID, FIndex, Fock, statistics, AbstractFockOperator, FOperator, BOperator, isnormalordered
export FCID, FockCoupling, fockcouplingnambus
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻
export @fc_str, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str
export Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb

"""
    ANNIHILATION

Indicate that the nambu index is ANNIHILATION.
"""
const ANNIHILATION = 1

"""
    CREATION

Indicate that the nambu index is CREATION.
"""
const CREATION = 2

"""
    MAJORANA

Indicate that the nambu index is MAJORANA.
"""
const MAJORANA = 0

"""
    FID <: IID

The Fock id.
"""
struct FID{ST} <: IID
    orbital::Int
    spin::Int
    nambu::Int
    function FID{ST}(orbital::Int, spin::Int, nambu::Int) where ST
        @assert ST∈(:f, :b) "FID error: wrong statistics."
        @assert nambu∈(0, 1, 2) "FID error: wrong input nambu($nambu)."
        new{ST}(orbital, spin, nambu)
    end
end
Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join(values(fid), ", ")
@inline Base.adjoint(fid::FID) = FID{statistics(fid)}(fid.orbital, fid.spin, fid.nambu==0 ? 0 : 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(fid::FID) = statistics(typeof(fid))
@inline statistics(::Type{FID{ST}}) where ST = ST

"""
    FID{ST}(; orbital::Int=1, spin::Int=1, nambu::Int=ANNIHILATION) where ST

Create a Fock id.
"""
@inline FID{ST}(; orbital::Int=1, spin::Int=1, nambu::Int=ANNIHILATION) where ST = FID{ST}(orbital, spin, nambu)

"""
    Fock{ST} <: Internal{FID{ST}}

The Fock interanl degrees of freedom.
"""
struct Fock{ST} <: Internal{FID{ST}}
    atom::Int
    norbital::Int
    nspin::Int
    nnambu::Int
    function Fock{ST}(atom::Int, norbital::Int, nspin::Int, nnambu::Int) where ST
        @assert ST∈(:f, :b) "Fock error: wrong statistics."
        @assert nnambu∈(1, 2) "Fock error: wrong input nnambu($nnambu)."
        new{ST}(atom, norbital, nspin, nnambu)
    end
end
@inline Base.Dims(fock::Fock) = (fock.norbital, fock.nspin, fock.nnambu)
@inline Base.CartesianIndex(fid::FID{ST}, fock::Fock{ST}) where ST = CartesianIndex(fid.orbital, fid.spin, fock.nnambu==1 ? 1 : fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2], fock.nnambu==1 ? 0 : index[3])
Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline statistics(fock::Fock) = statistics(typeof(fock))
@inline statistics(::Type{Fock{ST}}) where ST = ST

"""
    Fock{ST}(; atom::Int=1, norbital::Int=1, nspin::Int=2, nnambu::Int=2) where ST

Construct a Fock degrees of freedom.
"""
@inline Fock{ST}(; atom::Int=1, norbital::Int=1, nspin::Int=2, nnambu::Int=2) where ST = Fock{ST}(atom, norbital, nspin, nnambu)

"""
    FIndex{ST, S} <: Index{PID{S}, FID{ST}}

The Fock index.
"""
struct FIndex{ST, S} <: Index{PID{S}, FID{ST}}
    scope::S
    site::Int
    orbital::Int
    spin::Int
    nambu::Int
end
function FIndex{ST}(scope::S, site::Int, orbital::Int, spin::Int, nambu::Int) where {ST, S}
    @assert ST∈(:f, :b) "FIndex error: wrong statistics."
    @assert nambu∈(0, 1, 2) "FIndex error: wrong input nambu($nambu)."
    FIndex{ST, S}(scope, site, orbital, spin, nambu)
end
Base.show(io::IO, index::FIndex) = @printf io "FIndex{%s}(%s)" repr(statistics(index)) join(repr.(values(index)), ", ")
@inline Base.union(::Type{P}, ::Type{I}) where {P<:PID, I<:FID} = FIndex{statistics(I), fieldtype(P, :scope)}
@inline @generated function Base.replace(index::FIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){statistics(index)}($(exprs...)))
end
@inline statistics(index::FIndex) = statistics(typeof(index))
@inline statistics(::Type{<:FIndex{ST}}) where ST = ST

"""
    FIndex(pid::PID, fid::FID) -> FIndex

Construct a Fock index by a pid and an fid.
"""
@inline FIndex(pid::PID, fid::FID) = FIndex{statistics(fid)}(values(pid)..., values(fid)...)

"""
    script(::Val{:site}, index::FIndex; kwargs...) -> Int
    script(::Val{:orbital}, index::FIndex; kwargs...) -> Int
    script(::Val{:spinint}, index::FIndex; kwargs...) -> Int
    script(::Val{:spinsym}, index::FIndex; kwargs...) -> String
    script(::Val{:nambu}, index::FIndex; kwargs...) -> String

Get the required script of a Fock index.
"""
@inline script(::Val{:site}, index::FIndex; kwargs...) = index.site
@inline script(::Val{:orbital}, index::FIndex; kwargs...) = index.orbital
@inline script(::Val{:spinint}, index::FIndex; kwargs...) = index.spin
@inline script(::Val{:spinsym}, index::FIndex; kwargs...) = index.spin==1 ? "↓" : index.spin==2 ? "↑" : error("script error: wrong spin.")
@inline script(::Val{:nambu}, index::FIndex; kwargs...) = index.nambu==CREATION ? "\\dagger" : ""

"""
    fdefaultlatex

The default LaTeX format for a fermionic oid.
"""
const fdefaultlatex = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('c')
@inline latexname(::Type{<:FIndex{:f}}) = Symbol("FIndex{:f}")
@inline latexname(::Type{<:OID{<:FIndex{:f}}}) = Symbol("OID{FIndex{:f}}")
latexformat(FIndex{:f}, fdefaultlatex)
latexformat(OID{<:FIndex{:f}}, fdefaultlatex)

"""
    bdefaultlatex

The default LaTeX format for a bosonic oid.
"""
const bdefaultlatex = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('b')
@inline latexname(::Type{<:FIndex{:b}}) = Symbol("FIndex{:b}")
@inline latexname(::Type{<:OID{<:FIndex{:b}}}) = Symbol("OID{FIndex{:b}}")
latexformat(FIndex{:b}, bdefaultlatex)
latexformat(OID{<:FIndex{:b}}, bdefaultlatex)

"""
    angle(id::OID{<:FIndex}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float}) -> Complex{Float}

Get the twist phase corresponding to a Fock oid.
"""
function Base.angle(id::OID{<:FIndex}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase = length(vectors)==1 ? 2pi*dot(decompose(id.icoord, vectors[1]), values) :
            length(vectors)==2 ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            length(vectors)==3 ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return id.index.nambu==ANNIHILATION ? phase : id.index.nambu==CREATION ? -phase : error("angle error: not supported Fock index.")
end

"""
    usualfockindextotuple

Indicate that the choosed fields are `(:scope, :site, :orbital, :spin)` when converting a Fock index to tuple.
"""
const usualfockindextotuple = OIDToTuple(:scope, :site, :orbital, :spin)

"""
    nambufockindextotuple

Indicate that the choosed fields are `(:scope, :nambu, :site, :orbital, :spin)` when converting a Fock index to tuple.
"""
const nambufockindextotuple = OIDToTuple(:scope, :nambu, :site, :orbital, :spin)

"""
    AbstractFockOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex}}} <: Operator{V, I}

Abstract type for all Fock operators.
"""
abstract type AbstractFockOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex}}} <: Operator{V, I}end

"""
    isnormalordered(opt::AbstractFockOperator) -> Bool

Judge whether a Fock operator is normal ordered.
"""
function isnormalordered(opt::AbstractFockOperator)
    flag = true
    for i = 1:rank(opt)
        flag && (opt.id[i].index.nambu == ANNIHILATION) && (flag = false)
        flag || (opt.id[i].index.nambu == CREATION) && return false
    end
    return true
end

"""
    FOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex{:f}}}} <: AbstractFockOperator{V, I}

Fermionic Fock operator.
"""
struct FOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex{:f}}}} <: AbstractFockOperator{V, I}
    value::V
    id::I
    FOperator(value::Number) = new{typeof(value), Tuple{}}(value, ())
    FOperator(value::Number, id::Tuple{}) = new{typeof(value), Tuple{}}(value, ())
    FOperator(value::Number, id::ID{AbstractCompositeOID{<:FIndex{:f}}}) = new{typeof(value), typeof(id)}(value, id)
end

"""
    statistics(opt::FOperator) -> Symbol
    statistics(::Type{<:FOperator}) -> Symbol

Get the statistics of FOperator.
"""
@inline statistics(opt::FOperator) = statistics(typeof(opt))
@inline statistics(::Type{<:FOperator}) = :f

"""
    permute(::Type{<:FOperator}, id₁::OID{<:FIndex{:f}}, id₂::OID{<:FIndex{:f}}) -> Tuple{Vararg{FOperator}}

Permute two fermionic oid and get the result.
"""
function permute(::Type{<:FOperator}, id₁::OID{<:FIndex{:f}}, id₂::OID{<:FIndex{:f}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        return (FOperator(1), FOperator(-1, ID(id₂, id₁)))
    else
        return (FOperator(-1, ID(id₂, id₁)),)
    end
end

"""
    *(f1::FOperator, f2::FOperator) -> Union{Nothing, FOperator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(f1::FOperator, f2::FOperator)
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{Element, Element}, f1, f2)
end

"""
    BOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex{:b}}}} <: AbstractFockOperator{V, I}

Bosonic Fock operator.
"""
struct BOperator{V<:Number, I<:ID{AbstractCompositeOID{<:FIndex{:b}}}} <: AbstractFockOperator{V, I}
    value::V
    id::I
    BOperator(value::Number) = new{typeof(value), Tuple{}}(value, ())
    BOperator(value::Number, id::Tuple{}) = new{typeof(value), Tuple{}}(value, ())
    BOperator(value::Number, id::ID{AbstractCompositeOID{<:FIndex{:b}}}) = new{typeof(value), typeof(id)}(value, id)
end

"""
    statistics(opt::BOperator) -> Symbol
    statistics(::Type{<:BOperator}) -> Symbol

Get the statistics of BOperator.
"""
@inline statistics(opt::BOperator) = statistics(typeof(opt))
@inline statistics(::Type{<:BOperator}) = :b

"""
    permute(::Type{<:BOperator}, id₁::OID{<:FIndex{:b}}, id₂::OID{<:FIndex{:b}}) -> Tuple{Vararg{BOperator}}

Permute two bosonic oid and get the result.
"""
function permute(::Type{<:BOperator}, id₁::OID{<:FIndex{:b}}, id₂::OID{<:FIndex{:b}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        if id₁.index.nambu == CREATION
            return (BOperator(1), BOperator(1, ID(id₂, id₁)))
        else
            return (BOperator(-1), BOperator(1, ID(id₂, id₁)))
        end
    else
        return (BOperator(1, ID(id₂, id₁)),)
    end
end

"""
    FCID{A<:Tuple, N<:Tuple} <: SimpleID

The id of the atoms and nambus part of a Fock coupling.
"""
struct FCID{A<:Tuple, N<:Tuple} <: SimpleID
    atoms::A
    nambus::N
end

"""
    FockCoupling{V, A<:Tuple, N<:Tuple, O<:Subscripts, S<:Subscripts, I<:Tuple{FCID, SubID, SubID}} <: Coupling{V, I}

Fock coupling.
"""
struct FockCoupling{V, A<:Tuple, N<:Tuple, O<:Subscripts, S<:Subscripts, I<:Tuple{FCID, SubID, SubID}} <: Coupling{V, I}
    value::V
    atoms::A
    nambus::N
    orbitals::O
    spins::S
    function FockCoupling(value::Number, atoms::Tuple, nambus::Tuple, orbitals::Subscripts, spins::Subscripts)
        @assert length(atoms)==length(nambus)==length(orbitals)==length(spins) "FockCoupling error: dismatched atoms, nambus orbitals and spins."
        fcid, obid, spid = FCID(atoms, nambus), SubID(orbitals), SubID(spins)
        V, A, N, O, S, I = typeof(value), typeof(atoms), typeof(nambus), typeof(orbitals), typeof(spins),Tuple{typeof(fcid), typeof(obid), typeof(spid)}
        new{V, A, N, O, S, I}(value, atoms, nambus, orbitals, spins)
    end
end
@inline parameternames(::Type{<:FockCoupling}) = (:value, :atoms, :nambus, :orbitals, :spins, :id)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:atoms}, ::Type{A}) where {A<:Tuple} = !isconcretetype(A)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:nambus}, ::Type{N}) where {N<:Tuple} = !isconcretetype(N)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:orbitals}, ::Type{O}) where {O<:Subscripts} = !isconcretetype(O)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:spins}, ::Type{S}) where {S<:Subscripts} = !isconcretetype(S)
@inline contentnames(::Type{<:FockCoupling}) = (:value, :id, :orbitals, :spins)
@inline getcontent(fc::FockCoupling, ::Val{:id}) = ID(FCID(fc.atoms, fc.nambus), SubID(fc.orbitals), SubID(fc.spins))
@inline function FockCoupling(value::Number, id::Tuple{FCID, SubID, SubID}, orbitals::Subscripts, spins::Subscripts)
    FockCoupling(value, id[1].atoms, id[1].nambus, orbitals, spins)
end
@inline rank(::Type{<:FockCoupling{V, A} where V}) where A = fieldcount(A)

"""
    FockCoupling{N}(value::Number=1;
        atoms::Union{NTuple{N, Int}, Nothing}=nothing,
        nambus::Union{NTuple{N, Int}, Nothing}=nothing,
        orbitals::Union{NTuple{N, Int}, Subscripts, Nothing}=nothing,
        spins::Union{NTuple{N, Int}, Subscripts, Nothing}=nothing
        ) where N
"""
function FockCoupling{N}(value::Number=1;
        atoms::Union{NTuple{N, Int}, Nothing}=nothing,
        nambus::Union{NTuple{N, Int}, Nothing}=nothing,
        orbitals::Union{Int, NTuple{N, Int}, Subscripts}=N,
        spins::Union{Int, NTuple{N, Int}, Subscripts}=N
        ) where N
    isnothing(atoms) && (atoms = ntuple(i->wildcard, Val(N)))
    isnothing(nambus) && (nambus = ntuple(i->wildcard, Val(N)))
    isa(orbitals, Subscripts) || (orbitals = Subscripts(orbitals))
    isa(spins, Subscripts) || (spins = Subscripts(spins))
    return FockCoupling(value, atoms, nambus, orbitals, spins)
end

"""
    show(io::IO, fc::FockCoupling)

Show a Fock coupling.
"""
function Base.show(io::IO, fc::FockCoupling)
    @printf io "FockCoupling{%s}(value=%s" rank(fc) decimaltostr(fc.value)
    any(fc.atoms .≠ wildcard) && @printf io ", atoms=[%s]" join(fc.atoms, " ")
    any(fc.nambus .≠ wildcard) && @printf io ", nambus=[%s]" join(fc.nambus, " ")
    any(fc.orbitals .≠ wildcard) && @printf io ", orbitals=%s" string(fc.orbitals)
    any(fc.spins .≠ wildcard) && @printf io ", spins=%s" string(fc.spins)
    @printf io ")"
end

"""
    repr(fc::FockCoupling) -> String

Get the repr representation of a Fock coupling.
"""
function Base.repr(fc::FockCoupling)
    contents = String[]
    any(fc.atoms .≠ wildcard) && push!(contents, @sprintf "sl[%s]" join(fc.atoms, " "))
    any(fc.nambus .≠ wildcard) && push!(contents, @sprintf "ph[%s]" join(fc.nambus, " "))
    any(fc.orbitals .≠ wildcard) && push!(contents, @sprintf "ob%s" fc.orbitals)
    any(fc.spins .≠ wildcard) && push!(contents, @sprintf "sp%s" fc.spins)
    result = decimaltostr(fc.value)
    length(contents)==0 && (result = @sprintf "%s {%s}" result rank(fc))
    length(contents)>0 && (result = @sprintf "%s %s" result join(contents, " ⊗ "))
    return result
end

"""
    *(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the multiplication between two Fock couplings.
"""
@inline function Base.:*(fc₁::FockCoupling, fc₂::FockCoupling)
    return FockCoupling(fc₁.value*fc₂.value, (fc₁.atoms..., fc₂.atoms...), (fc₁.nambus..., fc₂.nambus...), fc₁.orbitals*fc₂.orbitals, fc₁.spins*fc₂.spins)
end

"""
    ⊗(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the direct product between two Fock couplings.
"""
function ⊗(fc₁::FockCoupling, fc₂::FockCoupling)
    @assert rank(fc₁)==rank(fc₂) "⊗ error: dismatched rank."
    wildcards = NTuple{rank(fc₁), Symbol}(wildcard for i = 1:rank(fc₁))
    atoms = fockcouplingchoicerule(fc₁.atoms, fc₂.atoms, wildcards, :atoms)
    nambus = fockcouplingchoicerule(fc₁.nambus, fc₂.nambus, wildcards, :nambus)
    orbitals = fockcouplingchoicerule(fc₁.orbitals, fc₂.orbitals, wildcards, :orbitals)
    spins = fockcouplingchoicerule(fc₁.spins, fc₂.spins, wildcards, :spins)
    return FockCoupling(fc₁.value*fc₂.value, atoms, nambus, orbitals, spins)
end
function fockcouplingchoicerule(v₁, v₂, wildcards::Tuple{Vararg{Symbol}}, field::Symbol)
    t₁, t₂ = all(v₁ .== wildcards), all(v₂ .== wildcards)
    @assert t₁||t₂ "fockcouplingchoicerule error: dismatched $field ($v₁ and $v₂)."
    return t₁ ? v₂ : v₁
end

"""
    ⋅(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the dot product of two rank-2 Fock couplings.

A rank-2 FockCoupling can be considered as a matrix acting on the sublattice, orbital, spin and nambu spaces.
Therefore, the dot product here is defined as the multiplication between such matrices.
"""
function ⋅(fc₁::FockCoupling, fc₂::FockCoupling)
    @assert rank(fc₁)==rank(fc₂)==2 "⋅ error: input fockcouplings must be rank-2."
    attrpairs = []
    value = fc₁.value * fc₂.value
    wildcards = (wildcard, wildcard)
    for attrname in (:atoms, :orbitals, :spins, :nambus)
        v₁ = (getproperty(fc₁, attrname)[1], getproperty(fc₁, attrname)[2])
        v₂ = (getproperty(fc₂, attrname)[1], getproperty(fc₂, attrname)[2])
        if isa(v₁, NTuple{2, Int}) && isa(v₂, NTuple{2, Int})
            value = value * delta(v₁[2], v₂[1])
            push!(attrpairs, attrname=>(v₁[1], v₂[2]))
        else
            @assert v₁==wildcards && v₂==wildcards "⋅ error: dismatched $attrname ($v₁ and $v₂)."
        end
    end
    return FockCoupling{2}(value; attrpairs...)
end

"""
    expand(fc::FockCoupling, points::NTuple{R, Point}, focks::NTuple{R, Fock}, info::Val) where R -> Union{FCExpand, Tuple{}}

Expand a Fock coupling with the given set of points and Fock degrees of freedom.
"""
function expand(fc::FockCoupling, points::NTuple{R, Point}, focks::NTuple{R, Fock}, info::Val) where R
    @assert rank(fc)==R "expand error: dismatched rank."
    for (i, atom) in enumerate(fc.atoms)
       isa(atom, Int) && atom≠focks[i].atom && return ()
    end
    nambus = fockcouplingnambus(fc, NTuple{rank(fc), Int}(focks[i].nnambu for i = 1:rank(fc)), info)
    obexpands = collect(expand(fc.orbitals, NTuple{rank(fc), Int}(focks[i].norbital for i = 1:rank(fc))))
    spexpands = collect(expand(fc.spins, NTuple{rank(fc), Int}(focks[i].nspin for i = 1:rank(fc))))
    return FCExpand{totalstatitics(focks)}(fc.value, points, obexpands, spexpands, nambus)
end
@generated totalstatitics(focks::NTuple{R, Fock}) where R = Tuple(statistics(fieldtype(focks, i)) for i = 1:R)
@generated function fockcouplingnambus(fc::FockCoupling, ranges::Tuple{Vararg{Int}}, ::Val)
    @assert rank(fc)==fieldcount(ranges) "fockcouplingnambus error: dismatched rank."
    exprs = [:(isa(fc.nambus[$i], Int) ?
            (ranges[$i]==1 && fc.nambus[$i]==0 || ranges[$i]==2 && 0<fc.nambus[$i]<=2 ? fc.nambus[$i] : error("fockcouplingnambus error: nambu out of range.")) :
            (ranges[$i]==1 ? 0 : ($i)%2==1 ? CREATION : ANNIHILATION)) for i = 1:rank(fc)]
    return Expr(:tuple, exprs...)
end
struct FCExpand{STS, V, N, D, S} <: CartesianVectorSpace{Tuple{V, ID{OID{FIndex, SVector{D, Float}}, N}}}
    value::V
    points::NTuple{N, Point{D, PID{S}}}
    obexpands::Vector{NTuple{N, Int}}
    spexpands::Vector{NTuple{N, Int}}
    nambus::NTuple{N, Int}
    function FCExpand{STS}(value, points::NTuple{N, Point{D, PID{S}}}, obexpands::Vector{NTuple{N, Int}}, spexpands::Vector{NTuple{N, Int}}, nambus::NTuple{N, Int}) where {STS, N, D, S}
        new{STS, typeof(value), N, D, S}(value, points, obexpands, spexpands, nambus)
    end
end
@inline @generated function Base.eltype(::Type{FCExpand{STS, V, N, D, S}}) where {STS, V, N, D, S}
    return Tuple{V, Tuple{[OID{FIndex{STS[i], S}, SVector{D, Float}} for i = 1:N]...}}
end
@inline Base.Dims(fce::FCExpand) = (length(fce.obexpands), length(fce.spexpands))
@generated function Tuple(index::CartesianIndex{2}, fce::FCExpand{STS, V, N}) where {STS, V, N}
    exprs = []
    for i = 1:N
        ST = QuoteNode(STS[i])
        push!(exprs, quote
            pid, rcoord, icoord = fce.points[$i].pid, fce.points[$i].rcoord, fce.points[$i].icoord
            fid = FID{$ST}(fce.obexpands[index[1]][$i], fce.spexpands[index[2]][$i], fce.nambus[$i])
            OID(FIndex(pid, fid), rcoord, icoord)
        end)
    end
    return Expr(:tuple, :(fce.value), Expr(:tuple, exprs...))
end

@inline σᵅname(mode) = mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : mode=="ph" ? :nambus : error("σᵅname error: wrong input mode.")
"""
    σ⁰(mode::String) -> Couplings

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁰(mode::String)
    attrname = σᵅname(mode)
    (attrval₁, attrval₂) = mode=="ph" ? ((ANNIHILATION, CREATION), (CREATION, ANNIHILATION)) : ((1, 1), (2, 2))
    FockCoupling{2}(1; attrname=>attrval₁) + FockCoupling{2}(1; attrname=>attrval₂)
end

"""
    σˣ(mode::String) -> Couplings

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σˣ(mode::String)
    attrname = σᵅname(mode)
    (attrval₁, attrval₂) = mode=="ph" ? ((ANNIHILATION, ANNIHILATION), (CREATION, CREATION)) : ((1, 2), (2, 1))
    FockCoupling{2}(1; attrname=>attrval₁) + FockCoupling{2}(1; attrname=>attrval₂)
end

"""
    σʸ(mode::String) -> Couplings

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σʸ(mode::String)
    attrname = σᵅname(mode)
    (attrval₁, attrval₂) = mode=="ph" ? ((ANNIHILATION, ANNIHILATION), (CREATION, CREATION)) : ((1, 2), (2, 1))
    FockCoupling{2}(1im; attrname=>attrval₁) + FockCoupling{2}(-1im; attrname=>attrval₂)
end

"""
    σᶻ(mode::String) -> Couplings

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σᶻ(mode::String)
    attrname = σᵅname(mode)
    (attrval₁, attrval₂) = mode=="ph" ? ((ANNIHILATION, CREATION), (CREATION, ANNIHILATION)) : ((1, 1), (2, 2))
    FockCoupling{2}(-1; attrname=>attrval₁) + FockCoupling{2}(1; attrname=>attrval₂)
end

"""
    σ⁺(mode::String) -> Couplings

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁺(mode::String)
    attrname = σᵅname(mode)
    attrval = mode=="ph" ? (CREATION, CREATION) : (2, 1)
    Couplings(FockCoupling{2}(1; attrname=>attrval))
end

"""
    σ⁻(mode::String) -> Couplings

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁻(mode::String)
    attrname = σᵅname(mode)
    attrval = mode=="ph" ? (ANNIHILATION, ANNIHILATION) : (1, 2)
    Couplings(FockCoupling{2}(1; attrname=>attrval))
end

"""
    fc"..." -> FockCoupling

Construct a FockCoupling from a literal string.
"""
macro fc_str(str)
    pos=findfirst(x->x∈('{', 's', 'o', 'p'), str)
    @assert 0<pos<length(str) "@fc_str error: wrong input pattern."
    coeff = eval(Meta.parse(str[1:pos-1]))
    ps = strip(replace(str[pos:end], " ⊗ "=>"⊗"))
    if ps[1]=='{' && ps[end]=='}'
        N = parse(Int, ps[2:end-1])
        return FockCoupling{N}(coeff)
    else
        attrpairs = []
        N = nothing
        for component in split(ps, '⊗')
            attrpair, rank = fccomponent(component)
            isnothing(N) && (N = rank)
            @assert N==rank "@fc_str error: dismatched ranks."
            push!(attrpairs, attrpair)
        end
        @assert length(attrpairs)>0 "@fc_str error: wrong input pattern."
        return Expr(:call, :(FockCoupling{$N}), Expr(:parameters, attrpairs...), coeff)
    end
end
function fccomponent(str::AbstractString)
    attrname = σᵅname(str[1:2])
    expr = Meta.parse(str[3:end])
    if attrname∈(:atoms, :nambus)
        @assert expr.head∈(:hcat, :vect) "fccomponent error: wrong input pattern for atoms and nambus."
        attrvalue = Tuple(expr.args)
        N = length(attrvalue)
    else
        @assert expr.head∈(:call, :hcat, :vcat, :vect) "fccomponent error: wrong input pattern for orbitals and spins."
        attrvalue = subscriptsexpr(expr)
        N = length(attrvalue.args[end].args[2])
    end
    return Expr(:kw, attrname, attrvalue), N
end

"""
    σ⁰"sp" -> Couplings
    σ⁰"ob" -> Couplings
    σ⁰"sl" -> Couplings
    σ⁰"ph" -> Couplings

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁰_str(mode::String) σ⁰(mode) end

"""
    σˣ"sp" -> Couplings
    σˣ"ob" -> Couplings
    σˣ"sl" -> Couplings
    σˣ"ph" -> Couplings

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σˣ_str(mode::String) σˣ(mode) end

"""
    σʸ"sp" -> Couplings
    σʸ"ob" -> Couplings
    σʸ"sl" -> Couplings
    σʸ"ph" -> Couplings

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σʸ_str(mode::String) σʸ(mode) end

"""
    σᶻ"sp" -> Couplings
    σᶻ"ob" -> Couplings
    σᶻ"sl" -> Couplings
    σᶻ"ph" -> Couplings

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σᶻ_str(mode::String) σᶻ(mode) end

"""
    σ⁺"sp" -> Couplings
    σ⁺"ob" -> Couplings
    σ⁺"sl" -> Couplings
    σ⁺"ph" -> Couplings

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁺_str(mode::String) σ⁺(mode) end

"""
    σ⁻"sp" -> Couplings
    σ⁻"ob" -> Couplings
    σ⁻"sl" -> Couplings
    σ⁻"ph" -> Couplings

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁻_str(mode::String) σ⁻(mode) end

"""
    Onsite(id::Symbol, value::Any;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Onsite term.

Type alias for `Term{:Onsite, 2, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const Onsite{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, 2, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value::Any;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? FockCoupling{2}() : couplings)
    Term{:Onsite, 2}(id, value, 0, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Onsite}) = :st
@inline isHermitian(::Type{<:Onsite}) = nothing

"""
    Hopping(id::Symbol, value::Any, bondkind::Int=1;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Hopping term.

Type alias for `Term{:Hopping, 2, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const Hopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, 2, id, V, Int, C, A, M}
@inline function Hopping(id::Symbol, value::Any, bondkind::Int=1;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? FockCoupling{2}() : couplings)
    @assert bondkind ≠ 0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    Term{:Hopping, 2}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hopping}) = :hp
@inline isHermitian(::Type{<:Hopping}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Hopping}) = (1, 2)

"""
    Pairing(id::Symbol, value::Any, bondkind::Int=0;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Pairing term.

Type alias for `Term{:Pairing, 2, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const Pairing{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, 2, id, V, Int, C, A, M}
@inline function Pairing(id::Symbol, value::Any, bondkind::Int=0;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? FockCoupling{2}() : couplings)
    Term{:Pairing, 2}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Pairing}) = :pr
@inline isHermitian(::Type{<:Pairing}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Pairing}) = (1, 2)
function fockcouplingnambus(fc::FockCoupling, ranges::NTuple{2, Int}, ::Val{:Pairing})
    @assert rank(fc)==2 "fockcouplingnambus error: dismatched rank."
    @assert ranges==(2, 2) "fockcouplingnambus error: ranges for Pairing terms must be (2, 2)."
    n₁ = isa(fc.nambus[1], Int) ? (0<fc.nambus[1]<=2 ? fc.nambus[1] : error("fockcouplingnambus error: nambu out of range.")) : ANNIHILATION
    n₂ = isa(fc.nambus[2], Int) ? (0<fc.nambus[2]<=2 ? fc.nambus[2] : error("fockcouplingnambus error: nambu out of range.")) : ANNIHILATION
    return (n₁, n₂)
end
function expand!(operators::Operators, term::Pairing, bond::AbstractBond, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing)
    argtypes = Tuple{Operators, Term, AbstractBond, Config, Bool}
    invoke(expand!, argtypes, operators, term, bond, config, half; table=table)
    isa(bond, Bond) && invoke(expand!, argtypes, operators, term, reverse(bond), config, half; table=table)
    return operators
end

"""
    Hubbard(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Hubbard term.

Type alias for `Term{:Hubbard, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, 4, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "Hubbard error: only real values are allowed."
    Term{:Hubbard, 4}(id, value, 0, couplings=fc"1 ph[2 1 2 1] ⊗ sp[2 2 1 1]", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hubbard}) = :hb
@inline isHermitian(::Type{<:Hubbard}) = true

"""
    InterOrbitalInterSpin(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, 4, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "InterOrbitalInterSpin error: only real values are allowed."
    Term{:InterOrbitalInterSpin, 4}(id, value, 0, couplings=fc"1 ph[2 1 2 1] ⊗ ob[α α β β](α < β) ⊗ sp[σ₁ σ₁ σ₂ σ₂](σ₁ ≠ σ₂)", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:InterOrbitalInterSpin}) = :nons
@inline isHermitian(::Type{<:InterOrbitalInterSpin}) = true

"""
    InterOrbitalIntraSpin(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, 4, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "InterOrbitalIntraSpin error: only real values are allowed."
    Term{:InterOrbitalIntraSpin, 4}(id, value, 0, couplings=fc"1 ph[2 1 2 1] ⊗ ob[α α β β](α < β)", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:InterOrbitalIntraSpin}) = :noes
@inline isHermitian(::Type{<:InterOrbitalIntraSpin}) = true

"""
    SpinFlip(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Spin-flip term.

Type alias for `Term{:SpinFlip, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, 4, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:SpinFlip, 4}(id, value, 0, couplings=fc"1 ph[2 2 1 1] ⊗ ob[α β α β](α < β) ⊗ sp[2 1 1 2]", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:SpinFlip}) = :sf
@inline isHermitian(::Type{<:SpinFlip}) = false

"""
    PairHopping(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Pair-hopping term.

Type alias for `Term{:PairHopping, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, 4, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value::Any; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:PairHopping, 4}(id, value, 0, couplings=fc"1 ph[2 2 1 1] ⊗ ob[α α β β](α < β) ⊗ sp[2 1 1 2]", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:PairHopping}) = :ph
@inline isHermitian(::Type{<:PairHopping}) = false

"""
    Coulomb(id::Symbol, value::Any, bondkind::Int=1;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Coulomb term.

Type alias for `Term{:Coulomb, 4, id, V, Int, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
"""
const Coulomb{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, 4, id, V, Int, C, A, M}
@inline function Coulomb(id::Symbol, value::Any, bondkind::Int=1;
        couplings::Union{Function, Coupling, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? FockCoupling{2}()*FockCoupling{2}() : couplings)
    @assert bondkind≠0 "Coulomb error: input bondkind cannot be 0. Use `Hubbard/InterOrbitalInterSpin/InterOrbitalIntraSpin/SpinFlip/PairHopping` instead."
    Term{:Coulomb, 4}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Coulomb}) = :cl
@inline isHermitian(::Type{<:Coulomb}) = nothing
@inline termfactor(id::ID{OID, 4}, ::Val{:Coulomb}) = id[2]'==id[1] && id[4]'==id[3] ? 2 : 1
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Coulomb}) = (1, 1, 2, 2)

"""
    otype(T::Type{<:Term}, C::Type{<:Config{<:Fock{:f}}}, B::Type{<:AbstractBond})
    otype(T::Type{<:Term}, C::Type{<:Config{<:Fock{:b}}}, B::Type{<:AbstractBond})

Get the operator type of a Fock term.
"""
@inline otype(T::Type{<:Term}, C::Type{<:Config{<:Fock{:f}}}, B::Type{<:AbstractBond}) = FOperator{valtype(T), ID{oidtype(valtype(C), eltype(B), kind(T)|>Val), rank(T)}}
@inline otype(T::Type{<:Term}, C::Type{<:Config{<:Fock{:b}}}, B::Type{<:AbstractBond}) = BOperator{valtype(T), ID{oidtype(valtype(C), eltype(B), kind(T)|>Val), rank(T)}}

end # module
