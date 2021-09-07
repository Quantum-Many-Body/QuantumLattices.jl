module FockPackage

using LinearAlgebra: dot
using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, AbstractBond, Bond, decompose
using ..DegreesOfFreedom: SimpleIID, CompositeIID, SimpleInternal, CompositeInternal, InternalIndex, IIDSpace, IIDConstrain, ConstrainID, Subscript, subscriptexpr
using ..DegreesOfFreedom: Index, LaTeX, OID, AbstractCompositeOID, latexformat, OIDToTuple, Operator, Operators, Hilbert, Table, wildcard, diagonal
using ..Terms: Coupling, Couplings, @couplings, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, delta, decimaltostr
using ...Prerequisites.Traits: rawtype
using ...Mathematics.AlgebraOverFields: ID, Element

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: couplingcenters, abbr, termfactor
import ...Interfaces: rank, ⊗, ⋅, expand, expand!, permute
import ...Mathematics.VectorSpaces: shape, ndimshape

export ANNIHILATION, CREATION, MAJORANA, fdefaultlatex, bdefaultlatex, usualfockindextotuple, nambufockindextotuple
export FID, Fock, statistics, isnormalordered
export FockCoupling, fockcouplingnambus
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
    FID{T, O<:Union{Int, Symbol}, S<:Union{Int, Symbol}, N<:Union{Int, Symbol}} <: SimpleIID

The Fock id.
"""
struct FID{T, O<:Union{Int, Symbol}, S<:Union{Int, Symbol}, N<:Union{Int, Symbol}} <: SimpleIID
    orbital::O
    spin::S
    nambu::N
    function FID{T}(orbital::Union{Int, Symbol}, spin::Union{Int, Symbol}, nambu::Union{Int, Symbol}) where T
        @assert T∈(:f, :b, wildcard) "FID error: wrong statistics."
        isa(nambu, Int) && @assert nambu∈(0, 1, 2) "FID error: wrong input nambu($nambu)."
        new{T, typeof(orbital), typeof(spin), typeof(nambu)}(orbital, spin, nambu)
    end
end
Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join(values(fid), ", ")
@inline Base.adjoint(fid::FID{T, <:Union{Int, Symbol}, <:Union{Int, Symbol}, Int}) where T = FID{T}(fid.orbital, fid.spin, fid.nambu==0 ? 0 : 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(fid::FID) = statistics(typeof(fid))
@inline statistics(::Type{<:FID{T}}) where T = T

"""
    FID{T}(; orbital::Union{Int, Symbol}=1, spin::Union{Int, Symbol}=1, nambu::Union{Int, Symbol}=ANNIHILATION) where T

Create a Fock id.
"""
@inline FID{T}(; orbital::Union{Int, Symbol}=1, spin::Union{Int, Symbol}=1, nambu::Union{Int, Symbol}=ANNIHILATION) where T = FID{T}(orbital, spin, nambu)

"""
    Fock{T} <: SimpleInternal{FID{T, Int, Int, Int}}

The Fock interanl degrees of freedom.
"""
struct Fock{T} <: SimpleInternal{FID{T, Int, Int, Int}}
    norbital::Int
    nspin::Int
    nnambu::Int
    function Fock{T}(norbital::Int, nspin::Int, nnambu::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        @assert nnambu∈(1, 2) "Fock error: wrong input nnambu($nnambu)."
        new{T}(norbital, nspin, nnambu)
    end
end
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, fock.nnambu==1 ? (0:0) : (1:2))
@inline ndimshape(::Type{<:Fock}) = 3
@inline Base.CartesianIndex(fid::FID{T}, fock::Fock{T}) where T = CartesianIndex(fid.orbital, fid.spin, fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2], index[3])
Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline statistics(fock::Fock) = statistics(typeof(fock))
@inline statistics(::Type{Fock{T}}) where T = T
function Base.show(io::IO, fock::Fock)
    @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
end
@generated function shape(iidspace::IIDSpace{I, V}) where {I<:CompositeIID{<:Tuple{Vararg{FID}}}, V<:CompositeInternal{<:Tuple{Vararg{Fock}}}}
    @assert rank(I)==rank(V) "shape error: dismatched composite iid and composite internal space."
    Kind = Val(kind(iidspace))
    Expr(:tuple, [:(shape(IIDSpace(iidspace.iid[$i], iidspace.internal[InternalIndex($i)], $Kind); order=$i)...) for i = 1:rank(I)]...)
end
@inline function shape(iidspace::IIDSpace{<:FID, <:Fock}; order)
    obrange = fidshape(kind(iidspace)|>Val, :orbital|>Val, iidspace.iid.orbital, iidspace.internal.norbital)
    sprange = fidshape(kind(iidspace)|>Val, :spin|>Val, iidspace.iid.spin, iidspace.internal.nspin)
    phrange = fidshape(kind(iidspace)|>Val, :nambu|>Val, iidspace.iid.nambu, iidspace.internal.nnambu; order=order)
    return (obrange, sprange, phrange)
end
@inline fidshape(::Val, ::Val, ::Symbol, n::Int) = 1:n
@inline fidshape(::Val, ::Val{field}, v::Int, n::Int) where field = ((@assert 0<v<n+1 "shape error: $field out of range."); v:v)
@inline fidshape(::Val, ::Val{:nambu}, ::Symbol, n::Int; order) = n==1 ? (MAJORANA:MAJORANA) : order%2==1 ? (CREATION:CREATION) : (ANNIHILATION:ANNIHILATION)
@inline fidshape(::Val, ::Val{:nambu}, v::Int, n::Int; order) = ((@assert n==1 ? v==0 : 0<v<n+1 "shape error: nambu out of range."); v:v)

"""
    Fock{T}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where T

Construct a Fock degrees of freedom.
"""
@inline Fock{T}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where T = Fock{T}(norbital, nspin, nnambu)

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
@inline latexname(::Type{<:Index{<:AbstractPID, <:FID{:f}}}) = Symbol("Index{AbstractPID, FID{:f}}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:f}}}")
latexformat(Index{<:AbstractPID, <:FID{:f}}, fdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}, fdefaultlatex)

"""
    bdefaultlatex

The default LaTeX format for a bosonic oid.
"""
const bdefaultlatex = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('b')
@inline latexname(::Type{<:Index{<:AbstractPID, <:FID{:b}}}) = Symbol("Index{AbstractPID, FID{:b}}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:b}}}")
latexformat(Index{<:AbstractPID, <:FID{:b}}, bdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}, bdefaultlatex)

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
    permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:f}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:f}}}) -> Tuple{Vararg{Operator}}

Permute two fermionic oid and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:f}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:f}}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        return (Operator(1), Operator(-1, ID(id₂, id₁)))
    else
        return (Operator(-1, ID(id₂, id₁)),)
    end
end

"""
    *(f1::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}}, f2::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}}) -> Union{Nothing, Operator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(f1::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}}, f2::Operator{<:Number, <:ID{AbstractCompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}})
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{Element, Element}, f1, f2)
end

"""
    permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:b}}}) -> Tuple{Vararg{Operator}}

Permute two bosonic oid and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:b}}})
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
    FockCoupling(value::Number,
        orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        spins::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        nambus::Union{NTuple{N, Int}, NTuple{N, Symbol}}
        ) where N

Fock coupling.

Type alias for "Coupling{V, I<:ID{FID}, C<:IIDConstrain, CI<:ConstrainID}."
"""
const FockCoupling{V, I<:ID{FID}, C<:IIDConstrain, CI<:ConstrainID} = Coupling{V, I, C, CI}
function FockCoupling(value::Number,
        orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        spins::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        nambus::Union{NTuple{N, Int}, NTuple{N, Symbol}}
        ) where N
    return Coupling(value, ID(FID{wildcard}, orbitals.pattern, spins.pattern, nambus), IIDConstrain((orbital=orbitals, spin=spins)))
end
function Base.show(io::IO, fc::FockCoupling)
    wildcards = ntuple(i->wildcard, Val(rank(fc)))
    @printf io "FockCoupling{%s}(value=%s" rank(fc) decimaltostr(fc.value)
    fc.cid.orbitals≠wildcards && @printf io ", orbitals=%s" repr(fc.constrain, 1:length(fc.constrain), :orbital)
    fc.cid.spins≠wildcards && @printf io ", spins=%s" repr(fc.constrain, 1:length(fc.constrain), :spin)
    fc.cid.nambus≠wildcards && @printf io ", nambus=[%s]" join(fc.cid.nambus, " ")
    @printf io ")"
end
function Base.repr(fc::FockCoupling)
    contents = String[]
    wildcards = ntuple(i->wildcard, Val(rank(fc)))
    fc.cid.orbitals≠wildcards && push!(contents, @sprintf "ob%s" repr(fc.constrain, 1:length(fc.constrain), :orbital))
    fc.cid.spins≠wildcards && push!(contents, @sprintf "sp%s" repr(fc.constrain, 1:length(fc.constrain), :spin))
    fc.cid.nambus≠wildcards && push!(contents, @sprintf "ph[%s]" join(fc.cid.nambus, " "))
    result = decimaltostr(fc.value)
    length(contents)==0 && (result = @sprintf "%s {%s}" result rank(fc))
    length(contents)>0 && (result = @sprintf "%s %s" result join(contents, " ⊗ "))
    return result
end

"""
    FockCoupling{N}(value::Number=1;
        orbitals::Union{NTuple{N, Int}, Subscript}=Subscript(N),
        spins::Union{NTuple{N, Int}, Subscript}=Subscript(N),
        nambus::Union{NTuple{N, Int}, NTuple{N, Symbol}}=ntuple(i->wildcard, Val(N))
        ) where N
"""
function FockCoupling{N}(value::Number=1;
        orbitals::Union{NTuple{N, Int}, Subscript}=Subscript(N),
        spins::Union{NTuple{N, Int}, Subscript}=Subscript(N),
        nambus::Union{NTuple{N, Int}, NTuple{N, Symbol}}=ntuple(i->wildcard, Val(N))
        ) where N
    isa(orbitals, Subscript) || (orbitals = Subscript(orbitals))
    isa(spins, Subscript) || (spins = Subscript(spins))
    return FockCoupling(value, orbitals, spins, nambus)
end

"""
    ⊗(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the direct product between two Fock couplings.
"""
function ⊗(fc₁::FockCoupling, fc₂::FockCoupling)
    @assert rank(fc₁)==rank(fc₂) "⊗ error: dismatched rank."
    @assert length(fc₁.constrain)==length(fc₂.constrain)==1 "⊗ error: not supported."
    wildcards = ntuple(i->wildcard, Val(rank(fc₁)))
    orbitals = fockcouplingchoicerule(first(fc₁.constrain).orbital, first(fc₂.constrain).orbital, wildcards, :orbitals)
    spins = fockcouplingchoicerule(first(fc₁.constrain).spin, first(fc₂.constrain).spin, wildcards, :spins)
    nambus = fockcouplingchoicerule(fc₁.cid.nambus, fc₂.cid.nambus, wildcards, :nambus)
    return FockCoupling(fc₁.value*fc₂.value, orbitals, spins, nambus)
end
function fockcouplingchoicerule(v₁, v₂, wildcards::Tuple{Vararg{Symbol}}, field::Symbol)
    t₁, t₂ = convert(Tuple, v₁)==wildcards, convert(Tuple, v₂)==wildcards
    @assert t₁||t₂ "fockcouplingchoicerule error: dismatched $field ($v₁ and $v₂)."
    return t₁ ? v₂ : v₁
end

"""
    ⋅(fc₁::FockCoupling, fc₂::FockCoupling) -> FockCoupling

Get the dot product of two rank-2 Fock couplings.

A rank-2 FockCoupling can be considered as a matrix acting on the orbital, spin and nambu spaces.
Therefore, the dot product here is defined as the multiplication between such matrices.
"""
function ⋅(fc₁::FockCoupling, fc₂::FockCoupling)
    @assert rank(fc₁)==rank(fc₂)==2 "⋅ error: input fockcouplings must be rank-2."
    value = fc₁.value * fc₂.value
    attrpairs = Pair{Symbol, NTuple{2, Int}}[]
    for attrname in (:orbitals, :spins, :nambus)
        v₁ = getproperty(fc₁.cid, attrname)
        v₂ = getproperty(fc₂.cid, attrname)
        if isa(v₁, NTuple{2, Int}) && isa(v₂, NTuple{2, Int})
            value = value * delta(v₁[2], v₂[1])
            push!(attrpairs, attrname=>(v₁[1], v₂[2]))
        else
            @assert v₁==(wildcard, wildcard) && v₂==(wildcard, wildcard) "⋅ error: dismatched $attrname ($v₁ and $v₂)."
        end
    end
    return FockCoupling{2}(value; attrpairs...)
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
macro fc_str(str::String)
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
        N = length(expr.args)
        attrvalue = Tuple(expr.args)
    else
        @assert expr.head∈(:call, :hcat, :vect) "fccomponent error: wrong input pattern for orbitals and spins."
        N = expr.head∈(:hcat, :vect) ? length(expr.args) : length(expr.args[1].args)
        attrvalue = subscriptexpr(expr)
    end
    return Expr(:kw, attrname, attrvalue), N
end

"""
    Onsite(id::Symbol, value;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Onsite term.

Type alias for `Term{:Onsite, 2, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Onsite{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, 2, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? @couplings(FockCoupling{2}()) : couplings)
    Term{:Onsite, 2}(id, value, 0, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Onsite}) = :st
@inline isHermitian(::Type{<:Onsite}) = nothing

"""
    Hopping(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Hopping term.

Type alias for `Term{:Hopping, 2, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, 2, id, V, Int, C, A, M}
@inline function Hopping(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? @couplings(FockCoupling{2}()) : couplings)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    Term{:Hopping, 2}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hopping}) = :hp
@inline isHermitian(::Type{<:Hopping}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Hopping}) = (1, 2)

"""
    Pairing(id::Symbol, value, bondkind::Int=0;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Pairing term.

Type alias for `Term{:Pairing, 2, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Pairing{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, 2, id, V, Int, C, A, M}
@inline function Pairing(id::Symbol, value, bondkind::Int=0;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? @couplings(FockCoupling{2}()) : couplings)
    Term{:Pairing, 2}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Pairing}) = :pr
@inline isHermitian(::Type{<:Pairing}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Pairing}) = (1, 2)
@inline fidshape(::Val{:Pairing}, ::Val{:nambu}, ::Symbol, n::Int; order) = ((@assert n==2 "range error: nnambu must be 2 for Pairing."); ANNIHILATION:ANNIHILATION)
function expand!(operators::Operators, term::Pairing, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    argtypes = Tuple{Operators, Term, AbstractBond, Hilbert, Bool}
    invoke(expand!, argtypes, operators, term, bond, hilbert, half; table=table)
    isa(bond, Bond) && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert, half; table=table)
    return operators
end

"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Hubbard term.

Type alias for `Term{:Hubbard, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, 4, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "Hubbard error: only real values are allowed."
    Term{:Hubbard, 4}(id, value, 0, couplings=@couplings(fc"1 sp[2 2 1 1] ⊗ ph[2 1 2 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hubbard}) = :hb
@inline isHermitian(::Type{<:Hubbard}) = true

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, 4, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "InterOrbitalInterSpin error: only real values are allowed."
    Term{:InterOrbitalInterSpin, 4}(id, value, 0, couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ sp[σ₁ σ₁ σ₂ σ₂](σ₁ ≠ σ₂) ⊗ ph[2 1 2 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:InterOrbitalInterSpin}) = :nons
@inline isHermitian(::Type{<:InterOrbitalInterSpin}) = true

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, 4, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert value==value' "InterOrbitalIntraSpin error: only real values are allowed."
    Term{:InterOrbitalIntraSpin, 4}(id, value, 0, couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ ph[2 1 2 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:InterOrbitalIntraSpin}) = :noes
@inline isHermitian(::Type{<:InterOrbitalIntraSpin}) = true

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Spin-flip term.

Type alias for `Term{:SpinFlip, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, 4, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:SpinFlip, 4}(id, value, 0, couplings=@couplings(fc"1 ob[α β α β](α < β) ⊗ sp[2 1 1 2] ⊗ ph[2 2 1 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:SpinFlip}) = :sf
@inline isHermitian(::Type{<:SpinFlip}) = false

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Pair-hopping term.

Type alias for `Term{:PairHopping, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, 4, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:PairHopping, 4}(id, value, 0, couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ sp[2 1 1 2] ⊗ ph[2 2 1 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:PairHopping}) = :ph
@inline isHermitian(::Type{<:PairHopping}) = false

"""
    Coulomb(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Coulomb term.

Type alias for `Term{:Coulomb, 4, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Coulomb{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, 4, id, V, Int, C, A, M}
@inline function Coulomb(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings, Nothing}=nothing,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    couplings = TermCouplings(isnothing(couplings) ? @couplings(FockCoupling{2}()*FockCoupling{2}()) : couplings)
    @assert bondkind≠0 "Coulomb error: input bondkind cannot be 0. Use `Hubbard/InterOrbitalInterSpin/InterOrbitalIntraSpin/SpinFlip/PairHopping` instead."
    Term{:Coulomb, 4}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Coulomb}) = :cl
@inline isHermitian(::Type{<:Coulomb}) = nothing
@inline termfactor(id::ID{OID, 4}, ::Val{:Coulomb}) = id[2]'==id[1] && id[4]'==id[3] ? 2 : 1
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Coulomb}) = (1, 1, 2, 2)

end # module
