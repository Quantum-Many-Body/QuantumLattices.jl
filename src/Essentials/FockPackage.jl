module FockPackage

using StaticArrays: SVector
using LinearAlgebra: dot
using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, AbstractBond, Point, Bond, decompose
using ..DegreesOfFreedom: SimpleIID, SimpleInternal, Index, LaTeX, OID, AbstractCompositeOID, latexformat, OIDToTuple, Operator, Operators, Hilbert, Table
using ..Terms: wildcard, Subscripts, SubID, Coupling, Couplings, couplingpoints, couplinginternals, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, delta, decimaltostr
using ...Prerequisites.Traits: rawtype
using ...Mathematics.AlgebraOverFields: SimpleID, ID, Element
using ...Mathematics.VectorSpaces: CartesianVectorSpace

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: subscriptsexpr, nonconstrain, couplingcenters, abbr, termfactor
import ...Interfaces: rank, ⊗, ⋅, expand, expand!, permute
import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent

export ANNIHILATION, CREATION, MAJORANA, fdefaultlatex, bdefaultlatex, usualfockindextotuple, nambufockindextotuple
export FID, Fock, statistics, isnormalordered
export FCID, FockCoupling, fockcouplingnambus
export @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str
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
    FID <: SimpleIID

The Fock id.
"""
struct FID{ST} <: SimpleIID
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
    Fock{ST} <: SimpleInternal{FID{ST}}

The Fock interanl degrees of freedom.
"""
struct Fock{ST} <: SimpleInternal{FID{ST}}
    norbital::Int
    nspin::Int
    nnambu::Int
    function Fock{ST}(norbital::Int, nspin::Int, nnambu::Int) where ST
        @assert ST∈(:f, :b) "Fock error: wrong statistics."
        @assert nnambu∈(1, 2) "Fock error: wrong input nnambu($nnambu)."
        new{ST}(norbital, nspin, nnambu)
    end
end
@inline Base.Dims(fock::Fock) = (fock.norbital, fock.nspin, fock.nnambu)
@inline Base.CartesianIndex(fid::FID{ST}, fock::Fock{ST}) where ST = CartesianIndex(fid.orbital, fid.spin, fock.nnambu==1 ? 1 : fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2], fock.nnambu==1 ? 0 : index[3])
Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline statistics(fock::Fock) = statistics(typeof(fock))
@inline statistics(::Type{Fock{ST}}) where ST = ST

"""
    Fock{ST}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where ST

Construct a Fock degrees of freedom.
"""
@inline Fock{ST}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where ST = Fock{ST}(norbital, nspin, nnambu)

"""
    script(::Val{:site}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:orbital}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:spinint}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:spinsym}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> String
    script(::Val{:nambu}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> String

Get the required script of a Fock index.
"""
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.pid.site
@inline script(::Val{:orbital}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.orbital
@inline script(::Val{:spinint}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.spin
@inline script(::Val{:spinsym}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.spin==1 ? "↓" : index.iid.spin==2 ? "↑" : error("script error: wrong spin.")
@inline script(::Val{:nambu}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.nambu==CREATION ? "\\dagger" : ""

"""
    fdefaultlatex

The default LaTeX format for a fermionic oid.
"""
const fdefaultlatex = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('c')
@inline latexname(::Type{<:Index{<:AbstractPID, FID{:f}}}) = Symbol("Index{AbstractPID, FID{:f}}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:f}}}")
latexformat(Index{<:AbstractPID, FID{:f}}, fdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}, fdefaultlatex)

"""
    bdefaultlatex

The default LaTeX format for a bosonic oid.
"""
const bdefaultlatex = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('b')
@inline latexname(::Type{<:Index{<:AbstractPID, FID{:b}}}) = Symbol("Index{AbstractPID, FID{:b}}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, FID{:b}}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:b}}}")
latexformat(Index{<:AbstractPID, FID{:b}}, bdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, FID{:b}}}, bdefaultlatex)

"""
    angle(id::OID{<:Index{<:AbstractPID, <:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Complex{<:Number}

Get the twist phase corresponding to a Fock oid.
"""
function Base.angle(id::OID{<:Index{<:AbstractPID, <:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    datatype = promote_type(eltype(values), Float)
    phase = length(vectors)==1 ? 2*convert(datatype, pi)*dot(decompose(id.icoord, vectors[1]), values) :
            length(vectors)==2 ? 2*convert(datatype, pi)*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            length(vectors)==3 ? 2*convert(datatype, pi)*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return id.index.iid.nambu==ANNIHILATION ? phase : id.index.iid.nambu==CREATION ? -phase : error("angle error: not supported Fock index.")
end

"""
    usualfockindextotuple

Indicate that the choosed fields are `(:site, :orbital, :spin)` when converting a Fock index to tuple.
"""
const usualfockindextotuple = OIDToTuple(:site, :orbital, :spin)

"""
    nambufockindextotuple

Indicate that the choosed fields are `(:nambu, :site, :orbital, :spin)` when converting a Fock index to tuple.
"""
const nambufockindextotuple = OIDToTuple(:nambu, :site, :orbital, :spin)

"""
    isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID}}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID}}}})
    flag = true
    for i = 1:rank(opt)
        flag && (opt.id[i].index.iid.nambu == ANNIHILATION) && (flag = false)
        flag || (opt.id[i].index.iid.nambu == CREATION) && return false
    end
    return true
end

"""
    permute(id₁::OID{<:Index{<:AbstractPID, FID{:f}}}, id₂::OID{<:Index{<:AbstractPID, FID{:f}}}) -> Tuple{Vararg{Operator}}

Permute two fermionic oid and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, FID{:f}}}, id₂::OID{<:Index{<:AbstractPID, FID{:f}}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        return (Operator(1), Operator(-1, ID(id₂, id₁)))
    else
        return (Operator(-1, ID(id₂, id₁)),)
    end
end

"""
    *(f1::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}}}, f2::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}}}) -> Union{Nothing, Operator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(f1::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}}}, f2::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, FID{:f}}}}})
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{Element, Element}, f1, f2)
end

"""
    permute(id₁::OID{<:Index{<:AbstractPID, FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, FID{:b}}}) -> Tuple{Vararg{Operator}}

Permute two bosonic oid and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, FID{:b}}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        if id₁.index.iid.nambu == CREATION
            return (Operator(1), Operator(1, ID(id₂, id₁)))
        else
            return (Operator(-1), Operator(1, ID(id₂, id₁)))
        end
    else
        return (Operator(1, ID(id₂, id₁)),)
    end
end

"""
    FCID{N<:Tuple} <: SimpleID

The id of the nambus part of a Fock coupling.
"""
struct FCID{N<:Tuple} <: SimpleID
    nambus::N
end

"""
    FockCoupling{V, N<:Tuple, O<:Subscripts, S<:Subscripts, I<:Tuple{FCID, SubID, SubID}} <: Coupling{V, I}

Fock coupling.
"""
struct FockCoupling{V, N<:Tuple, O<:Subscripts, S<:Subscripts, I<:Tuple{FCID, SubID, SubID}} <: Coupling{V, I}
    value::V
    nambus::N
    orbitals::O
    spins::S
    function FockCoupling(value::Number, nambus::Tuple, orbitals::Subscripts, spins::Subscripts)
        @assert length(nambus)==length(orbitals)==length(spins) "FockCoupling error: dismatched nambus, orbitals and spins."
        fcid, obid, spid = FCID(nambus), SubID(orbitals), SubID(spins)
        V, N, O, S, I = typeof(value), typeof(nambus), typeof(orbitals), typeof(spins),Tuple{typeof(fcid), typeof(obid), typeof(spid)}
        new{V, N, O, S, I}(value, nambus, orbitals, spins)
    end
end
@inline parameternames(::Type{<:FockCoupling}) = (:value, :nambus, :orbitals, :spins, :id)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:nambus}, ::Type{N}) where {N<:Tuple} = !isconcretetype(N)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:orbitals}, ::Type{O}) where {O<:Subscripts} = !isconcretetype(O)
@inline isparameterbound(::Type{<:FockCoupling}, ::Val{:spins}, ::Type{S}) where {S<:Subscripts} = !isconcretetype(S)
@inline contentnames(::Type{<:FockCoupling}) = (:value, :id, :orbitals, :spins)
@inline getcontent(fc::FockCoupling, ::Val{:id}) = ID(FCID(fc.nambus), SubID(fc.orbitals), SubID(fc.spins))
@inline function FockCoupling(value::Number, id::Tuple{FCID, SubID, SubID}, orbitals::Subscripts, spins::Subscripts)
    FockCoupling(value, id[1].nambus, orbitals, spins)
end
@inline rank(::Type{<:FockCoupling{V, N} where V}) where N = fieldcount(N)

"""
    FockCoupling{N}(value::Number=1;
        nambus::NTuple{N, Union{Int, Symbol}}=ntuple(i->wildcard, Val(N)),
        orbitals::Union{NTuple{N, Int}, Subscripts}=Subscripts(N),
        spins::Union{NTuple{N, Int}, Subscripts}=Subscripts(N)
        ) where N
"""
function FockCoupling{N}(value::Number=1;
        nambus::NTuple{N, Union{Int, Symbol}}=ntuple(i->wildcard, Val(N)),
        orbitals::Union{NTuple{N, Int}, Subscripts}=Subscripts(N),
        spins::Union{NTuple{N, Int}, Subscripts}=Subscripts(N)
        ) where N
    isa(orbitals, Subscripts) || (orbitals = Subscripts(orbitals))
    isa(spins, Subscripts) || (spins = Subscripts(spins))
    return FockCoupling(value, nambus, orbitals, spins)
end

"""
    show(io::IO, fc::FockCoupling)

Show a Fock coupling.
"""
function Base.show(io::IO, fc::FockCoupling)
    @printf io "FockCoupling{%s}(value=%s" rank(fc) decimaltostr(fc.value)
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
    return FockCoupling(fc₁.value*fc₂.value, (fc₁.nambus..., fc₂.nambus...), fc₁.orbitals*fc₂.orbitals, fc₁.spins*fc₂.spins)
end

"""
    ⊗(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the direct product between two Fock couplings.
"""
function ⊗(fc₁::FockCoupling, fc₂::FockCoupling)
    @assert rank(fc₁)==rank(fc₂) "⊗ error: dismatched rank."
    wildcards = NTuple{rank(fc₁), Symbol}(wildcard for i = 1:rank(fc₁))
    nambus = fockcouplingchoicerule(fc₁.nambus, fc₂.nambus, wildcards, :nambus)
    orbitals = fockcouplingchoicerule(fc₁.orbitals, fc₂.orbitals, wildcards, :orbitals)
    spins = fockcouplingchoicerule(fc₁.spins, fc₂.spins, wildcards, :spins)
    return FockCoupling(fc₁.value*fc₂.value, nambus, orbitals, spins)
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
    for attrname in (:orbitals, :spins, :nambus)
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
    expand(fc::FockCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val) -> FCExpand

Expand a Fock coupling with the given bond and the Hilbert space.
"""
function expand(fc::FockCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)
    points = couplingpoints(fc, bond, info)
    focks = couplinginternals(fc, bond, hilbert, info)
    @assert rank(fc)==length(points)==length(focks) "expand error: dismatched rank."
    nambus = fockcouplingnambus(fc, NTuple{rank(fc), Int}(focks[i].nnambu for i = 1:rank(fc)), info)
    obexpands = collect(expand(fc.orbitals, NTuple{rank(fc), Int}(focks[i].norbital for i = 1:rank(fc))))
    spexpands = collect(expand(fc.spins, NTuple{rank(fc), Int}(focks[i].nspin for i = 1:rank(fc))))
    return FCExpand{statistics(focks)}(fc.value, points, obexpands, spexpands, nambus)
end
@generated statistics(focks::NTuple{R, Fock}) where R = Tuple(statistics(fieldtype(focks, i)) for i = 1:R)
@generated function fockcouplingnambus(fc::FockCoupling, ranges::Tuple{Vararg{Int}}, ::Val)
    @assert rank(fc)==fieldcount(ranges) "fockcouplingnambus error: dismatched rank."
    exprs = [:(isa(fc.nambus[$i], Int) ?
            (ranges[$i]==1 && fc.nambus[$i]==0 || ranges[$i]==2 && 0<fc.nambus[$i]<=2 ? fc.nambus[$i] : error("fockcouplingnambus error: nambu out of range.")) :
            (ranges[$i]==1 ? 0 : ($i)%2==1 ? CREATION : ANNIHILATION)) for i = 1:rank(fc)]
    return Expr(:tuple, exprs...)
end
struct FCExpand{STS, V, N, D, P<:AbstractPID, DT<:Number} <: CartesianVectorSpace{Tuple{V, ID{OID{<:Index, SVector{D, DT}}, N}}}
    value::V
    points::NTuple{N, Point{D, P, DT}}
    obexpands::Vector{NTuple{N, Int}}
    spexpands::Vector{NTuple{N, Int}}
    nambus::NTuple{N, Int}
    function FCExpand{STS}(value, points::NTuple{N, Point{D, P, DT}}, obexpands::Vector{NTuple{N, Int}}, spexpands::Vector{NTuple{N, Int}}, nambus::NTuple{N, Int}) where {STS, N, D, P<:AbstractPID, DT<:Number}
        new{STS, typeof(value), N, D, P, DT}(value, points, obexpands, spexpands, nambus)
    end
end
@inline @generated function Base.eltype(::Type{FCExpand{STS, V, N, D, P, DT}}) where {STS, V, N, D, P<:AbstractPID, DT<:Number}
    return Tuple{V, Tuple{[OID{Index{P, FID{STS[i]}}, SVector{D, DT}} for i = 1:N]...}}
end
@inline Base.Dims(fce::FCExpand) = (length(fce.obexpands), length(fce.spexpands))
@generated function Tuple(index::CartesianIndex{2}, fce::FCExpand{STS, V, N}) where {STS, V, N}
    exprs = []
    for i = 1:N
        ST = QuoteNode(STS[i])
        push!(exprs, quote
            pid, rcoord, icoord = fce.points[$i].pid, fce.points[$i].rcoord, fce.points[$i].icoord
            fid = FID{$ST}(fce.obexpands[index[1]][$i], fce.spexpands[index[2]][$i], fce.nambus[$i])
            OID(Index(pid, fid), rcoord, icoord)
        end)
    end
    return Expr(:tuple, :(fce.value), Expr(:tuple, exprs...))
end

"""
    σ⁰"sp" -> Couplings
    σ⁰"ob" -> Couplings
    σ⁰"ph" -> Couplings

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σ⁰_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σ⁰_str error: wrong input mode."
    mode=="sp" && return FockCoupling{2}(1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    mode=="ob" && return FockCoupling{2}(1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    return FockCoupling{2}(1, nambus=(ANNIHILATION, CREATION)) + FockCoupling{2}(1, nambus=(CREATION, ANNIHILATION))
end

"""
    σˣ"sp" -> Couplings
    σˣ"ob" -> Couplings
    σˣ"ph" -> Couplings

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σˣ_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σˣ_str error: wrong input mode."
    mode=="sp" && return FockCoupling{2}(1, spins=(1, 2)) + FockCoupling{2}(1, spins=(2, 1))
    mode=="ob" && return FockCoupling{2}(1, orbitals=(1, 2)) + FockCoupling{2}(1, orbitals=(2, 1))
    return FockCoupling{2}(1, nambus=(ANNIHILATION, ANNIHILATION)) + FockCoupling{2}(1, nambus=(CREATION, CREATION))
end

"""
    σʸ"sp" -> Couplings
    σʸ"ob" -> Couplings
    σʸ"ph" -> Couplings

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σʸ_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σʸ_str error: wrong input mode."
    mode=="sp" && return FockCoupling{2}(1im, spins=(1, 2)) + FockCoupling{2}(-1im, spins=(2, 1))
    mode=="ob" && return FockCoupling{2}(1im, orbitals=(1, 2)) + FockCoupling{2}(-1im, orbitals=(2, 1))
    return FockCoupling{2}(1im, nambus=(ANNIHILATION, ANNIHILATION)) + FockCoupling{2}(-1im, nambus=(CREATION, CREATION))
end

"""
    σᶻ"sp" -> Couplings
    σᶻ"ob" -> Couplings
    σᶻ"ph" -> Couplings

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σᶻ_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σᶻ_str error: wrong input mode."
    mode=="sp" && return FockCoupling{2}(-1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    mode=="ob" && return FockCoupling{2}(-1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    return FockCoupling{2}(-1, nambus=(ANNIHILATION, CREATION)) + FockCoupling{2}(1, nambus=(CREATION, ANNIHILATION))
end

"""
    σ⁺"sp" -> Couplings
    σ⁺"ob" -> Couplings
    σ⁺"ph" -> Couplings

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σ⁺_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σ⁺_str error: wrong input mode."
    mode=="sp" && return Couplings(FockCoupling{2}(1, spins=(2, 1)))
    mode=="ob" && return Couplings(FockCoupling{2}(1, orbitals=(2, 1)))
    return Couplings(FockCoupling{2}(1, nambus=(CREATION, CREATION)))
end

"""
    σ⁻"sp" -> Couplings
    σ⁻"ob" -> Couplings
    σ⁻"ph" -> Couplings

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob") or particle-holes("ph").
"""
macro σ⁻_str(mode::String)
    @assert mode∈("sp", "ob", "ph") "@σ⁻_str error: wrong input mode."
    mode=="sp" && return Couplings(FockCoupling{2}(1, spins=(1, 2)))
    mode=="ob" && return Couplings(FockCoupling{2}(1, orbitals=(1, 2)))
    return Couplings(FockCoupling{2}(1, nambus=(ANNIHILATION, ANNIHILATION)))
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
@inline σᵅname(mode) = mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="ph" ? :nambus : error("σᵅname error: wrong input mode.")
function fccomponent(str::AbstractString)
    attrname = σᵅname(str[1:2])
    expr = Meta.parse(str[3:end])
    if attrname==:nambus
        @assert expr.head∈(:hcat, :vect) "fccomponent error: wrong input pattern for nambus."
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
function expand!(operators::Operators, term::Pairing, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    argtypes = Tuple{Operators, Term, AbstractBond, Hilbert, Bool}
    invoke(expand!, argtypes, operators, term, bond, hilbert, half; table=table)
    isa(bond, Bond) && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert, half; table=table)
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

end # module
