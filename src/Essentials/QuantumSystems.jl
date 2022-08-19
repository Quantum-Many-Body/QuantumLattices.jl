module QuantumSystems

using LinearAlgebra: dot, norm
using StaticArrays: SVector
using Printf: @printf, @sprintf
using ..QuantumOperators: ID, OperatorProd, Operator, Operators, LaTeX, latexformat
using ..Spatials: AbstractPID, AbstractBond, Point, Bond, rcoord, pidtype
using ..DegreesOfFreedom: SimpleIID, CompositeIID, SimpleInternal, CompositeInternal
using ..DegreesOfFreedom: Index, CompositeOID, OID, OIDToTuple, Hilbert, Table
using ..DegreesOfFreedom: Subscripts, Subscript, IIDSpace, subscriptexpr, wildcard, diagonal
using ..DegreesOfFreedom: Coupling, Couplings, @couplings, couplinginternals, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind, dtype
using ...Interfaces: decompose, dimension
using ...Prerequisites: Float, atol, rtol, delta, decimaltostr, concatenate
using ...Prerequisites.Traits: efficientoperations, rawtype, getcontent
using ...Prerequisites.VectorSpaces: VectorSpace, VectorSpaceStyle, VectorSpaceCartesian

import ..QuantumOperators: script, latexname, ishermitian, optype
import ..DegreesOfFreedom: statistics, couplingcenters, abbr, termfactor
import ...Interfaces: rank, ⊗, ⋅, expand, expand!, permute
import ...Prerequisites.VectorSpaces: shape

export majorana, annihilation, creation, flatex, blatex
export FID, Fock, FockCoupling, Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb, FockTerm
export isnormalordered, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str

export slatex
export SID, Spin, SpinCoupling, SpinTerm
export totalspin, @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

export nlatex
export NID, Phonon, PhononCoupling, PhononKinetic, PhononPotential, PhononTerm
export @kinetic_str, @potential_str

# Canonical fermionic/bosonic systems and hardcore bosonic systems
"""
    majorana

Indicate that the nambu index is majorana.
"""
const majorana = 0

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
@inline Base.:(==)(fid₁::FID, fid₂::FID) = statistics(fid₁)==statistics(fid₂) && ==(efficientoperations, fid₁, fid₂)
@inline Base.isequal(fid₁::FID, fid₂::FID) = isequal(statistics(fid₁), statistics(fid₂)) && isequal(efficientoperations, fid₁, fid₂)
Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join(ntuple(i->getfield(fid, i), Val(fieldcount(typeof(fid)))), ", ")
@inline Base.adjoint(fid::FID{T, <:Union{Int, Symbol}, <:Union{Int, Symbol}, Int}) where T = FID{T}(fid.orbital, fid.spin, fid.nambu==0 ? 0 : 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(::Type{<:FID{T}}) where T = T
@inline FID(iid::FID, ::CompositeInternal) = iid

"""
    FID{T}(; orbital::Union{Int, Symbol}=1, spin::Union{Int, Symbol}=1, nambu::Union{Int, Symbol}=annihilation) where T

Create a Fock id.
"""
@inline FID{T}(; orbital::Union{Int, Symbol}=1, spin::Union{Int, Symbol}=1, nambu::Union{Int, Symbol}=annihilation) where T = FID{T}(orbital, spin, nambu)

"""
    Fock{T} <: SimpleInternal{FID{T, Int, Int, Int}}

The Fock internal degrees of freedom.
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
@inline Base.eltype(::Type{Fock}) = (FID{T, Int, Int, Int} where T)
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, fock.nnambu==1 ? (0:0) : (1:2))
@inline Base.CartesianIndex(fid::FID{T}, fock::Fock{T}) where T = CartesianIndex(fid.orbital, fid.spin, fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2], index[3])
Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
function Base.show(io::IO, fock::Fock)
    @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
end
@inline Base.match(::Type{<:FID{wildcard}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false

@inline function CompositeIID(iidspace::IIDSpace{I, V}) where {I<:CompositeIID{<:Tuple{Vararg{FID}}}, V<:CompositeInternal{:⊗, <:Tuple{Vararg{Fock}}}}
    @assert rank(I)==rank(V) "CompositeIID error: mismatched composite iid and composite internal space."
    return CompositeIID(map(
        (order, fid, fock)->FID{statistics(fid)}(fid.orbital, fid.spin, fidnambu(Val(kind(iidspace)), fid.nambu, fock.nnambu, order)),
        ntuple(i->i, Val(rank(I))),
        iidspace.iid.contents,
        iidspace.internal.contents
        ))
end
@inline fidnambu(::Val, ::Symbol, nnambu::Int, order::Int) = nnambu==1 ? majorana : order%2==1 ? creation : annihilation
@inline fidnambu(::Val, nambu::Int, nnambu::Int, ::Int) = ((@assert nnambu==1 ? nambu==0 : 0<nambu<nambu+1 "fidnambu error: nambu out of range."); nambu)

@inline function shape(iidspace::IIDSpace{<:FID, <:Fock})
    @assert isa(iidspace.iid.nambu, Int) "shape error: nambu not determined."
    obrange = fidshape(iidspace.iid.orbital, iidspace.internal.norbital)
    sprange = fidshape(iidspace.iid.spin, iidspace.internal.nspin)
    phrange = iidspace.iid.nambu:iidspace.iid.nambu
    return (obrange, sprange, phrange)
end
@inline fidshape(::Symbol, n::Int) = 1:n
@inline fidshape(v::Int, n::Int) where field = ((@assert 0<v<n+1 "shape error: $field out of range."); v:v)

"""
    Fock{T}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where T

Construct a Fock degrees of freedom.
"""
@inline Fock{T}(; norbital::Int=1, nspin::Int=2, nnambu::Int=2) where T = Fock{T}(norbital, nspin, nnambu)

"""
    script(::Val{:site}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:orbital}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:spint}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> Int
    script(::Val{:spsym}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> String
    script(::Val{:nambu}, index::Index{<:AbstractPID, <:FID}; kwargs...) -> String

Get the required script of a Fock index.
"""
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.pid.site
@inline script(::Val{:orbital}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.orbital
@inline script(::Val{:spint}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.spin
@inline script(::Val{:spsym}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.spin==1 ? "↓" : index.iid.spin==2 ? "↑" : error("script error: wrong spin.")
@inline script(::Val{:nambu}, index::Index{<:AbstractPID, <:FID}; kwargs...) = index.iid.nambu==creation ? "\\dagger" : ""

"""
    flatex

The default LaTeX format for a fermionic oid.
"""
const flatex = LaTeX{(:nambu,), (:site, :orbital, :spsym)}('c')
@inline latexname(::Type{<:Index{<:AbstractPID, <:FID{:f}}}) = Symbol("Index{AbstractPID, FID{:f}}")
@inline latexname(::Type{<:CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}) = Symbol("CompositeOID{Index{AbstractPID, FID{:f}}}")
latexformat(Index{<:AbstractPID, <:FID{:f}}, flatex)
latexformat(CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}, flatex)

"""
    blatex

The default LaTeX format for a bosonic oid.
"""
const blatex = LaTeX{(:nambu,), (:site, :orbital, :spsym)}('b')
@inline latexname(::Type{<:Index{<:AbstractPID, <:FID{:b}}}) = Symbol("Index{AbstractPID, FID{:b}}")
@inline latexname(::Type{<:CompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}}) = Symbol("CompositeOID{Index{AbstractPID, FID{:b}}}")
latexformat(Index{<:AbstractPID, <:FID{:b}}, blatex)
latexformat(CompositeOID{<:Index{<:AbstractPID, <:FID{:b}}}, blatex)

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
    return id.index.iid.nambu==annihilation ? phase : id.index.iid.nambu==creation ? -phase : error("angle error: not supported Fock index.")
end

"""
    isnormalordered(opt::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID}}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID}}}})
    flag = true
    for i = 1:rank(opt)
        flag && (opt.id[i].index.iid.nambu == annihilation) && (flag = false)
        flag || (opt.id[i].index.iid.nambu == creation) && return false
    end
    return true
end

"""
    permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:f}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:f}}}) -> Tuple{Vararg{Operator}}

Permute two fermionic oids and get the result.
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
    *(  f1::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}}
        ) -> Union{Nothing, Operator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(
        f1::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{CompositeOID{<:Index{<:AbstractPID, <:FID{:f}}}}}
        )
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{OperatorProd, OperatorProd}, f1, f2)
end

"""
    permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:b}}}) -> Tuple{Vararg{Operator}}

Permute two bosonic oids and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, <:FID{:b}}}, id₂::OID{<:Index{<:AbstractPID, <:FID{:b}}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        if id₁.index.iid.nambu == creation
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

Type alias for `Coupling{V, I<:ID{FID}, C<:Subscripts}`.
"""
const FockCoupling{V, I<:ID{FID}, C<:Subscripts} = Coupling{V, I, C}
function FockCoupling(value::Number,
        orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        spins::Subscript{<:NTuple{N, Union{Int, Symbol}}},
        nambus::Union{NTuple{N, Int}, NTuple{N, Symbol}}
        ) where N
    return Coupling(value, ID(FID{wildcard}, orbitals.pattern, spins.pattern, nambus), Subscripts((orbital=orbitals, spin=spins)))
end
function Base.show(io::IO, fc::FockCoupling)
    wildcards = ntuple(i->wildcard, Val(rank(fc)))
    @printf io "FockCoupling{%s}(value=%s" rank(fc) decimaltostr(fc.value)
    fc.cid.orbitals≠wildcards && @printf io ", orbitals=%s" repr(fc.subscripts, 1:length(fc.subscripts), :orbital)
    fc.cid.spins≠wildcards && @printf io ", spins=%s" repr(fc.subscripts, 1:length(fc.subscripts), :spin)
    fc.cid.nambus≠wildcards && @printf io ", nambus=[%s]" join(fc.cid.nambus, " ")
    @printf io ")"
end
function Base.repr(fc::FockCoupling)
    contents = String[]
    wildcards = ntuple(i->wildcard, Val(rank(fc)))
    fc.cid.orbitals≠wildcards && push!(contents, @sprintf "ob%s" repr(fc.subscripts, 1:length(fc.subscripts), :orbital))
    fc.cid.spins≠wildcards && push!(contents, @sprintf "sp%s" repr(fc.subscripts, 1:length(fc.subscripts), :spin))
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

Construct a Fock coupling.
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
    @assert rank(fc₁)==rank(fc₂) "⊗ error: mismatched rank."
    @assert length(fc₁.subscripts)==length(fc₂.subscripts)==1 "⊗ error: not supported."
    wildcards = ntuple(i->wildcard, Val(rank(fc₁)))
    orbitals = fockcouplingchoicerule(first(fc₁.subscripts).orbital, first(fc₂.subscripts).orbital, wildcards, :orbitals)
    spins = fockcouplingchoicerule(first(fc₁.subscripts).spin, first(fc₂.subscripts).spin, wildcards, :spins)
    nambus = fockcouplingchoicerule(fc₁.cid.nambus, fc₂.cid.nambus, wildcards, :nambus)
    return FockCoupling(fc₁.value*fc₂.value, orbitals, spins, nambus)
end
function fockcouplingchoicerule(v₁, v₂, wildcards::Tuple{Vararg{Symbol}}, field::Symbol)
    t₁, t₂ = Tuple(v₁)==wildcards, Tuple(v₂)==wildcards
    @assert t₁||t₂ "fockcouplingchoicerule error: mismatched $field ($v₁ and $v₂)."
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
            @assert v₁==(wildcard, wildcard) && v₂==(wildcard, wildcard) "⋅ error: mismatched $attrname ($v₁ and $v₂)."
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
    return FockCoupling{2}(1, nambus=(annihilation, creation)) + FockCoupling{2}(1, nambus=(creation, annihilation))
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
    return FockCoupling{2}(1, nambus=(annihilation, annihilation)) + FockCoupling{2}(1, nambus=(creation, creation))
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
    return FockCoupling{2}(1im, nambus=(annihilation, annihilation)) + FockCoupling{2}(-1im, nambus=(creation, creation))
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
    return FockCoupling{2}(-1, nambus=(annihilation, creation)) + FockCoupling{2}(1, nambus=(creation, annihilation))
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
    return Couplings(FockCoupling{2}(1, nambus=(creation, creation)))
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
    return Couplings(FockCoupling{2}(1, nambus=(annihilation, annihilation)))
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
            @assert N==rank "@fc_str error: mismatched ranks."
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
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Onsite{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    Term{:Onsite}(id, value, 0, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Onsite}) = :st
@inline ishermitian(::Type{<:Onsite}) = nothing

"""
    Hopping(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Hopping term.

Type alias for `Term{:Hopping, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, id, V, Int, C, A, M}
@inline function Hopping(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    Term{:Hopping}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hopping}) = :hp
@inline ishermitian(::Type{<:Hopping}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Hopping}) = (1, 2)

"""
    Pairing(id::Symbol, value, bondkind::Int=0;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Pairing term.

Type alias for `Term{:Pairing, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Pairing{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, id, V, Int, C, A, M}
@inline function Pairing(id::Symbol, value, bondkind::Int=0;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    Term{:Pairing}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Pairing}) = :pr
@inline ishermitian(::Type{<:Pairing}) = false
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Pairing}) = (1, 2)
@inline fidnambu(::Val{:Pairing}, ::Symbol, nnambu::Int, order::Int) = ((@assert nnambu==2 "fidnambu error: nnambu must be 2 for Pairing."); annihilation)
function expand!(operators::Operators, term::Pairing, bond::AbstractBond, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    argtypes = Tuple{Operators, Term, AbstractBond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half, table=table)
    isa(bond, Bond) && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half, table=table)
    return operators
end

"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Hubbard term.

Type alias for `Term{:Hubbard, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:Hubbard}(id, value, 0, couplings=@couplings(fc"1 sp[2 2 1 1] ⊗ ph[2 1 2 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Hubbard}) = :hb
@inline ishermitian(::Type{<:Hubbard}) = true

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:InterOrbitalInterSpin}(id, value, 0;
        couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ sp[σ₁ σ₁ σ₂ σ₂](σ₁ ≠ σ₂) ⊗ ph[2 1 2 1]"),
        amplitude=amplitude,
        modulate=modulate
        )
end
@inline abbr(::Type{<:InterOrbitalInterSpin}) = :nons
@inline ishermitian(::Type{<:InterOrbitalInterSpin}) = true

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:InterOrbitalIntraSpin}(id, value, 0, couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ ph[2 1 2 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:InterOrbitalIntraSpin}) = :noes
@inline ishermitian(::Type{<:InterOrbitalIntraSpin}) = true

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:SpinFlip}(id, value, 0, couplings=@couplings(fc"1 ob[α β α β](α < β) ⊗ sp[2 1 1 2] ⊗ ph[2 2 1 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:SpinFlip}) = :sf
@inline ishermitian(::Type{<:SpinFlip}) = false

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    Term{:PairHopping}(id, value, 0, couplings=@couplings(fc"1 ob[α α β β](α < β) ⊗ sp[2 1 1 2] ⊗ ph[2 2 1 1]"), amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:PairHopping}) = :ph
@inline ishermitian(::Type{<:PairHopping}) = false

"""
    Coulomb(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()*FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Coulomb term.

Type alias for `Term{:Coulomb, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const Coulomb{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, id, V, Int, C, A, M}
@inline function Coulomb(id::Symbol, value, bondkind::Int=1;
        couplings::Union{Function, Couplings}=@couplings(FockCoupling{2}()*FockCoupling{2}()),
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    @assert bondkind≠0 "Coulomb error: input bondkind cannot be 0. Use `Hubbard/InterOrbitalInterSpin/InterOrbitalIntraSpin/SpinFlip/PairHopping` instead."
    Term{:Coulomb}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:Coulomb}) = :cl
@inline ishermitian(::Type{<:Coulomb}) = nothing
@inline termfactor(id::ID{OID, 4}, ::Val{:Coulomb}) = id[2]'==id[1] && id[4]'==id[3] ? 2 : 1
@inline couplingcenters(::FockCoupling, ::Bond, ::Val{:Coulomb}) = (1, 1, 2, 2)

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

"""
    SID{S, O<:Union{Int, Symbol}} <: SimpleIID

The spin id.
"""
struct SID{S, O<:Union{Int, Symbol}, T<:Union{Char, Symbol}} <: SimpleIID
    orbital::O
    tag::T
    function SID{S}(orbital::Union{Int, Symbol}, tag::Union{Char, Symbol}) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) || S==wildcard "SID error: not supported spin($S)."
        isa(tag, Char) && @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
        new{S, typeof(orbital), typeof(tag)}(orbital, tag)
    end
end
@inline SID{S}(tag::Char) where S = SID{S}(1, tag)
@inline Base.:(==)(sid₁::SID, sid₂::SID) = totalspin(sid₁)==totalspin(sid₂) && ==(efficientoperations, sid₁, sid₂)
@inline Base.isequal(sid₁::SID, sid₂::SID) = isequal(totalspin(sid₁), totalspin(sid₂)) && isequal(efficientoperations, sid₁, sid₂)
@inline Base.adjoint(sid::SID) = SID{totalspin(sid)}(sid.orbital, sidajointmap[sid.tag])
Base.show(io::IO, sid::SID) = @printf io "SID{%s}(%s)" totalspin(sid) join(map(repr, ntuple(i->getfield(sid, i), Val(fieldcount(typeof(sid))))), ", ")
@inline @generated function Base.replace(sid::SID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(sid, $name))) for name in QuoteNode.(fieldnames(sid))]
    return :(rawtype(typeof(sid)){totalspin(sid)}($(exprs...)))
end
@inline totalspin(sid::SID) = totalspin(typeof(sid))
@inline totalspin(::Type{<:SID{S}}) where S = S
@inline statistics(::Type{<:SID}) = :b
@inline SID(iid::SID, ::CompositeInternal) = iid

"""
    SID{S}(tag::Union{Char, Symbol}; orbital::Union{Int, Symbol}=1) where S

Create a spin id.
"""
@inline SID{S}(tag::Union{Char, Symbol}; orbital::Union{Int, Symbol}=1) where S = SID{S}(orbital, tag)

"""
    Matrix(sid::SID{S, <:Union{Int, Symbol}, Char}, dtype::Type{<:Number}=Complex{Float}) where S -> Matrix{dtype}

Get the matrix representation of a sid.
"""
function Base.Matrix(sid::SID{S, <:Union{Int, Symbol}, Char}, dtype::Type{<:Number}=Complex{Float}) where S
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

"""
    Spin{S} <: SimpleInternal{SID{S, Int, Char}}

The spin internal degrees of freedom.
"""
struct Spin{S} <: SimpleInternal{SID{S, Int, Char}}
    norbital::Int
    function Spin{S}(norbital::Int) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
        new{S}(norbital)
    end
end
@inline Base.eltype(::Type{Spin}) = (SID{S, Int, Char} where S)
@inline shape(sp::Spin) = (1:sp.norbital, 1:length(sidtagmap))
@inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sid.orbital, sidseqmap[sid.tag])
@inline SID(index::CartesianIndex{2}, sp::Spin) = SID{totalspin(sp)}(index[1], sidtagmap[index[2]])
Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S
function Base.show(io::IO, spin::Spin)
    @printf io "%s{%s}(%s)" spin|>typeof|>nameof totalspin(spin) join(("$name=$(getfield(spin, name))" for name in spin|>typeof|>fieldnames), ", ")
end
@inline Base.match(::Type{<:SID{wildcard}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SID{S}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SID{S₁}}, ::Type{<:Spin{S₂}}) where {S₁, S₂} = false
@inline sidrange(::Symbol) = 1:5
@inline sidrange(::Symbol, n::Int) = 1:n
@inline sidrange(tag::Char) = (pos=sidseqmap[tag]; pos:pos)
@inline sidrange(v::Int, n::Int) = ((@assert 0<v<n+1 "shape error: orbital out of range."); v:v)
@inline shape(iidspace::IIDSpace{<:SID, <:Spin}) = (sidrange(iidspace.iid.orbital, iidspace.internal.norbital), sidrange(iidspace.iid.tag))

"""
    Spin{S}(; norbital::Int=1) where S

Construct a spin degrees of freedom.
"""
@inline Spin{S}(; norbital::Int=1) where S = Spin{S}(norbital)

"""
    script(::Val{:site}, index::Index{<:AbstractPID, <:SID}; kwargs...) -> Int
    script(::Val{:orbital}, index::Index{<:AbstractPID, <:SID}; kwargs...) -> Int
    script(::Val{:tag}, index::Index{<:AbstractPID, <:SID}; kwargs...) -> Char

Get the required script of a spin oid.
"""
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:SID}; kwargs...) = index.pid.site
@inline script(::Val{:orbital}, index::Index{<:AbstractPID, <:SID}; kwargs...) = index.iid.orbital
@inline script(::Val{:tag}, index::Index{<:AbstractPID, <:SID}; kwargs...) = index.iid.tag

"""
    slatex

The default LaTeX format for a spin oid.
"""
const slatex = LaTeX{(:tag,), (:site,)}('S')
@inline latexname(::Type{<:Index{<:AbstractPID, <:SID}}) = Symbol("Index{AbstractPID, SID}")
@inline latexname(::Type{<:CompositeOID{<:Index{<:AbstractPID, <:SID}}}) = Symbol("CompositeOID{Index{AbstractPID, SID}}")
latexformat(Index{<:AbstractPID, <:SID}, slatex)
latexformat(CompositeOID{<:Index{<:AbstractPID, <:SID}}, slatex)

"""
    permute(id₁::OID{<:Index{<:AbstractPID, SID{S, Int, Char}}}, id₂::OID{<:Index{<:AbstractPID, SID{S, Int, Char}}}) where S -> Tuple{Vararg{Operator}}

Permute two spin oids and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, SID{S, Int, Char}}}, id₂::OID{<:Index{<:AbstractPID, SID{S, Int, Char}}}) where S
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if id₁.index.pid==id₂.index.pid && id₁.index.iid.orbital==id₂.index.iid.orbital && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        if id₁.index.iid.tag == 'x'
            id₂.index.iid.tag=='y' && return (Operator(+1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='z' && return (Operator(-1im, ID(permutesoid(id₁, 'y'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='+' && return (Operator(-1, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='-' && return (Operator(+1, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
        elseif id₁.index.iid.tag == 'y'
            id₂.index.iid.tag=='x' && return (Operator(-1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='z' && return (Operator(+1im, ID(permutesoid(id₁, 'x'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='+' && return (Operator(-1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='-' && return (Operator(-1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
        elseif id₁.index.iid.tag == 'z'
            id₂.index.iid.tag=='x' && return (Operator(+1im, ID(permutesoid(id₁, 'y'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='y' && return (Operator(-1im, ID(permutesoid(id₁, 'x'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='+' && return (Operator(+1, ID(id₂)), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='-' && return (Operator(-1, ID(id₂)), Operator(1, ID(id₂, id₁)))
        elseif id₁.index.iid.tag == '+'
            id₂.index.iid.tag=='x' && return (Operator(+1, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='y' && return (Operator(+1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='z' && return (Operator(-1, ID(id₁)), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='-' && return (Operator(+2, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
        elseif id₁.index.iid.tag == '-'
            id₂.index.iid.tag=='x' && return (Operator(-1, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='y' && return (Operator(1im, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='z' && return (Operator(+1, ID(id₁)), Operator(1, ID(id₂, id₁)))
            id₂.index.iid.tag=='+' && return (Operator(-2, ID(permutesoid(id₁, 'z'))), Operator(1, ID(id₂, id₁)))
        end
    else
        return (Operator(1, ID(id₂, id₁)),)
    end
end
@inline permutesoid(id::OID{<:Index{<:AbstractPID, <:SID}}, tag::Char) = replace(id, index=replace(id.index, iid=replace(id.index.iid, tag=tag)))

"""
    SpinCoupling(value::Number, tags::NTuple{N, Char}, orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}}) where N

Spin coupling.

Type alias for `Coupling{V, I<:ID{SID}, C<:Subscripts}`.
"""
const SpinCoupling{V, I<:ID{SID}, C<:Subscripts} = Coupling{V, I, C}
@inline function SpinCoupling(value::Number, tags::NTuple{N, Char}, orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}}) where N
    return Coupling(value, ID(SID{wildcard}, orbitals.pattern, tags), Subscripts((orbital=orbitals,)))
end
function Base.show(io::IO, sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    @printf io ", tags=%s" join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.cid.tags), "")
    sc.cid.orbitals≠ntuple(i->wildcard, Val(rank(sc))) && @printf io ", orbitals=%s" repr(sc.subscripts, 1:length(sc.subscripts), :orbital)
    @printf io ")"
end
function Base.repr(sc::SpinCoupling)
    result = [@sprintf "%s %s" decimaltostr(sc.value) join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.cid.tags), "")]
    sc.cid.orbitals≠ntuple(i->wildcard, Val(rank(sc))) && push!(result, @sprintf "ob%s" repr(sc.subscripts, 1:length(sc.subscripts), :orbital))
    return join(result, " ")
end

"""
    SpinCoupling(value::Number, tags::NTuple{N, Char}; orbitals::Union{NTuple{N, Int}, Subscript}=Subscript(N)) where N

Construct a spin coupling.
"""
function SpinCoupling(value::Number, tags::NTuple{N, Char}; orbitals::Union{NTuple{N, Int}, Subscript}=Subscript(N)) where N
    isa(orbitals, Subscript) || (orbitals = Subscript(orbitals))
    return SpinCoupling(value, tags, orbitals)
end

"""
    sc"..." -> SpinCoupling

Construct a SpinCoupling from a literal string.
"""
macro sc_str(str::String)
    fpos = findfirst(r"S[⁺⁻ˣʸᶻ⁰]", str).start
    lpos = 1 + lastindex(str) - findfirst(r"[⁺⁻ˣʸᶻ⁰]S", str|>reverse).start
    coeff = eval(Meta.parse(str[firstindex(str):prevind(str, fpos)]))
    tags = Tuple(sidreprevmap[tag] for tag in replace(str[thisind(str, fpos):thisind(str, lpos)], "S"=>""))
    orbitals = scorbitals(strip(str[nextind(str, lpos):end]), length(tags))
    return Expr(:call, :SpinCoupling, orbitals, coeff, tags)
end
function scorbitals(str::AbstractString, n::Int)
    length(str)==0 && return Expr(:kw, :orbitals, Subscript(n))
    @assert str[1:2]=="ob" "scorbitals error: wrong input pattern."
    expr = Meta.parse(str[3:end])
    @assert expr.head∈(:call, :hcat, :vect) "scorbitals error: wrong input pattern for orbitals."
    attrvalue = subscriptexpr(expr)
    return Expr(:kw, :orbitals, attrvalue)
end

"""
    heisenberg"ob[o₁ o₂]" -> Couplings
    heisenberg"xyz ob[o₁ o₂]" -> Couplings
    heisenberg"+-z ob[o₁ o₂]" -> Couplings

The Heisenberg couplings.
"""
macro heisenberg_str(str::String)
    slice = findfirst(r"xyz|\+\-z", str)
    mode = isnothing(slice) ? "+-z" : str[slice]
    isnothing(slice) || (str = strip(str[nextind(str, slice.stop):end]))
    @assert mode=="+-z" || mode=="xyz" "@heisenberg_str error: not supported mode($mode)."
    orbitals = scorbitals(str, 2)
    if mode == "+-z"
        return Expr(:call, :heisenbergpmz, orbitals)
    else
        return Expr(:call, :heisenbergxyz, orbitals)
    end
end
function heisenbergxyz(;orbitals=Subscript(2))
    return Couplings(
        SpinCoupling(1, ('x', 'x'), orbitals=orbitals),
        SpinCoupling(1, ('y', 'y'), orbitals=orbitals),
        SpinCoupling(1, ('z', 'z'), orbitals=orbitals)
    )
end
function heisenbergpmz(;orbitals=Subscript(2))
    return Couplings(
        SpinCoupling(1//2, ('+', '-'), orbitals=orbitals),
        SpinCoupling(1//2, ('-', '+'), orbitals=orbitals),
        SpinCoupling(1//1, ('z', 'z'), orbitals=orbitals)
    )
end

"""
    ising"x ob[o₁ o₂]" -> Couplings
    ising"y ob[o₁ o₂])" -> Couplings
    ising"z ob[o₁ o₂]" -> Couplings

The Ising couplings.
"""
macro ising_str(str::String)
    @assert str[1] ∈ ('x', 'y', 'z') "@ising_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@ising_str error: wrong input pattern."
    orbitals = scorbitals(str[3:end], 2)
    return :(Couplings(SpinCoupling(1, ($str[1], $str[1]), $orbitals)))
end

"""
    gamma"x ob[o₁ o₂]" -> Couplings
    gamma"y ob[o₁ o₂]" -> Couplings
    gamma"z ob[o₁ o₂]" -> Couplings

The Gamma couplings.
"""
macro gamma_str(str::String)
    @assert str[1] in ('x', 'y', 'z') "@gamma_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@gamma_str error: wrong input pattern."
    t₁, t₂ = str[1]=='x' ? ('y', 'z') : str[1]=='y' ? ('z', 'x') : ('x', 'y')
    orbitals = scorbitals(str[3:end], 2)
    return Expr(:call, :gamma, orbitals, t₁, t₂)
end
function gamma(t₁::Char, t₂::Char; orbitals=Subscript(2))
    return Couplings(
        SpinCoupling(1, (t₁, t₂), orbitals=orbitals),
        SpinCoupling(1, (t₂, t₁), orbitals=orbitals)
        )
end

"""
    dm"x ob[o₁ o₂]" -> Couplings
    dm"y ob[o₁ o₂]" -> Couplings
    dm"z ob[o₁ o₂]" -> Couplings

The DM couplings.
"""
macro dm_str(str::String)
    @assert str[1] in ('x', 'y', 'z') "@dm_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@dm_str error: wrong input pattern."
    orbitals = scorbitals(str[3:end], 2)
    t₁, t₂ = str[1]=='x' ? ('y', 'z') : str[1]=='y' ? ('z', 'x') : ('x', 'y')
    return Expr(:call, :dm, orbitals, t₁, t₂)
end
function dm(t₁::Char, t₂::Char; orbitals=Subscript(2))
    return Couplings(
        SpinCoupling(+1, (t₁, t₂), orbitals=orbitals),
        SpinCoupling(-1, (t₂, t₁), orbitals=orbitals)
        )
end

"""
    sˣ"ob[o]" -> Couplings
    sʸ"ob[o]" -> Couplings
    sᶻ"ob[o]" -> Couplings

The single Sˣ/Sʸ/Sᶻ coupling.
"""
macro sˣ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, scorbitals(str, 1), 1, ('x',))) end
macro sʸ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, scorbitals(str, 1), 1, ('y',))) end
macro sᶻ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, scorbitals(str, 1), 1, ('z',))) end

"""
    SpinTerm(id::Symbol, value::Any, bondkind::Any, couplings::Union{Function, Coupling, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        )

Spin term.

Type alias for `Term{:SpinTerm, id, V, B<:Any, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinTerm{id, V, B<:Any, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:SpinTerm, id, V, B, C, A, M}
@inline function SpinTerm(id::Symbol, value::Any, bondkind::Any, couplings::Union{Function, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    Term{:SpinTerm}(id, value, bondkind, couplings=couplings, amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:SpinTerm}) = :sp
@inline ishermitian(::Type{<:SpinTerm}) = true
@inline function couplingcenters(sc::SpinCoupling, ::Bond, ::Val{:SpinTerm})
    @assert rank(sc)%2==0 "couplingcenters error: the rank of the input spin coupling should be even."
    return ntuple(i->2-i%2, Val(rank(sc)))
end

# Phononic systems
"""
    NID{D<:Union{Char, Symbol}} <: SimpleIID

The phonon id.
"""
struct NID{D<:Union{Char, Symbol}} <: SimpleIID
    tag::Char
    dir::D
    function NID(tag::Char, dir::Union{Char, Symbol}=wildcard)
        @assert tag∈('p', 'u') "NID error: wrong tag($tag)."
        isa(dir, Char) && @assert dir∈('x', 'y', 'z') "NID error: wrong direction($dir)."
        new{typeof(dir)}(tag, dir)
    end
end
@inline Base.adjoint(pnid::NID) = pnid
@inline statistics(::Type{<:NID}) = :b
@inline NID(iid::NID, ::CompositeInternal) = iid

"""
    Phonon <: SimpleInternal{NID{Char}}

The phonon internal degrees of freedom.
"""
struct Phonon <: SimpleInternal{NID{Char}}
    ndir::Int
    function Phonon(ndir::Integer)
        @assert ndir∈(1, 2, 3) "Phonon error: wrong number of directions."
        new(ndir)
    end
end
@inline shape(pn::Phonon) = (1:2, 1:pn.ndir)
@inline Base.CartesianIndex(pnid::NID{Char}, ::Phonon) = CartesianIndex(pnid.tag=='u' ? 1 : 2, Int(pnid.dir)-Int('x')+1)
@inline NID(index::CartesianIndex{2}, ::Phonon) = NID(index[1]==1 ? 'u' : 'p', Char(Int('x')+index[2]-1))
@inline shape(iidspace::IIDSpace{NID{Symbol}, Phonon}) = (iidspace.iid.tag=='u' ? (1:1) : (2:2), 1:iidspace.internal.ndir)
@inline function shape(iidspace::IIDSpace{NID{Char}, Phonon})
    dir = Int(iidspace.iid.dir)-Int('x')+1
    @assert 0<dir<iidspace.internal.ndir+1 "shape error: dir out of range."
    return (iidspace.iid.tag=='u' ? (1:1) : (2:2), dir:dir)
end

"""
    script(::Val{:BD}, index::Index{<:AbstractPID, <:NID}, l::LaTeX) -> Char
    script(::Val{:BD}, oid::CompositeOID{<:Index{<:AbstractPID, <:NID}}, l::LaTeX) -> Char
    script(::Val{:site}, index::Index{<:AbstractPID, <:NID}; kwargs...) -> Int
    script(::Val{:dir}, index::Index{<:AbstractPID, <:NID}; kwargs...) -> Char

Get the required script of a phonon oid.
"""
@inline script(::Val{:BD}, index::Index{<:AbstractPID, <:NID}, l::LaTeX) = l.body[index.iid.tag]
@inline script(::Val{:BD}, oid::CompositeOID{<:Index{<:AbstractPID, <:NID}}, l::LaTeX) = l.body[getcontent(oid, :index).iid.tag]
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:NID}; kwargs...) = index.pid.site
@inline script(::Val{:dir}, index::Index{<:AbstractPID, <:NID}; kwargs...) = index.iid.dir

"""
    nlatex

The default LaTeX format for a phonon oid.
"""
const nlatex = LaTeX{(:dir,), (:site,)}(Dict('p'=>'p', 'u'=>'u'), "", "")
@inline latexname(::Type{<:Index{<:AbstractPID, <:NID}}) = Symbol("Index{AbstractPID, NID}")
@inline latexname(::Type{<:CompositeOID{<:Index{<:AbstractPID, <:NID}}}) = Symbol("CompositeOID{Index{AbstractPID, NID}}")
latexformat(Index{<:AbstractPID, <:NID}, nlatex)
latexformat(CompositeOID{<:Index{<:AbstractPID, <:NID}}, nlatex)

"""
    permute(id₁::OID{<:Index{<:AbstractPID, NID{Char}}}, id₂::OID{<:Index{<:AbstractPID, NID{Char}}}) -> Tuple{Vararg{Operator}}

Permute two phonon oids and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, NID{Char}}}, id₂::OID{<:Index{<:AbstractPID, NID{Char}}})
    if id₁.index.iid.dir==id₂.index.iid.dir && id₁.index.iid.tag≠id₂.index.iid.tag && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        if id₁.index.iid.tag=='u'
            return (Operator(1im), Operator(1, ID(id₂, id₁)))
        else
            return (Operator(-1im), Operator(1, ID(id₂, id₁)))
        end
    else
        return (Operator(1, ID(id₂, id₁)),)
    end
end

"""
    PhononCoupling(value::Number, tags::NTuple{N, Char}, dirs::Subscript{<:Union{NTuple{N, Char}, NTuple{N, Symbol}}}) where N

Phonon coupling.

Type alias for `Coupling{V<:Number, I<:ID{NID}, C<:Subscripts}`.
"""
const PhononCoupling{V<:Number, I<:ID{NID}, C<:Subscripts} = Coupling{V, I, C}
@inline function PhononCoupling(value::Number, tags::NTuple{N, Char}, dirs::Subscript{<:Union{NTuple{N, Char}, NTuple{N, Symbol}}}) where N
    return Coupling(value, ID(NID, tags, dirs.pattern), Subscripts((dir=dirs,)))
end
function Base.show(io::IO, pnc::PhononCoupling)
    @printf io "PhononCoupling(value=%s" decimaltostr(pnc.value)
    @printf io ", tags=[%s]" join(pnc.cid.tags, " ")
    pnc.cid.dirs≠ntuple(i->wildcard, Val(rank(pnc))) && @printf io ", dirs=%s" repr(pnc.subscripts, 1:length(pnc.subscripts), :dir)
    @printf io ")"
end
function Base.repr(pnc::PhononCoupling)
    result = [@sprintf "%s [%s]" decimaltostr(pnc.value) join(pnc.cid.tags, " ")]
    pnc.cid.dirs≠ntuple(i->wildcard, Val(rank(pnc))) && push!(result, @sprintf "dr%s" repr(pnc.subscripts, 1:length(pnc.subscripts), :dir))
    return join(result, " ")
end

"""
    PhononCoupling(value::Number, tags::NTuple{N, Char};
        dirs::Union{NTuple{N, Char}, NTuple{N, Symbol}, Subscript}=Subscript(N)
        ) where N

Construct a phonon coupling.
"""
function PhononCoupling(value::Number, tags::NTuple{N, Char};
        dirs::Union{NTuple{N, Char}, NTuple{N, Symbol}, Subscript}=Subscript(N)
        ) where N
    isa(dirs, Subscript) || (dirs = Subscript(dirs))
    return PhononCoupling(value, tags, dirs)
end

"""
    expand(pnc::PhononCoupling{<:Number, <:ID{NID{Symbol}}}, bond::Bond, hilbert::Hilbert, info::Val{:PhononPotential}) -> PPExpand

Expand the default phonon potential coupling on a given bond.
"""
function expand(pnc::PhononCoupling{<:Number, <:ID{NID{Symbol}}}, bond::Bond, hilbert::Hilbert, info::Val{:PhononPotential})
    R̂ = rcoord(bond)/norm(rcoord(bond))
    @assert pnc.cid.tags==('u', 'u') "expand error: wrong tags of phonon coupling."
    @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of phonon coupling."
    pn₁, pn₂ = couplinginternals(pnc, bond, hilbert, info)
    @assert pn₁.ndir==pn₂.ndir==length(R̂) "expand error: mismatched number of directions."
    return PPExpand(R̂, (bond.epoint, bond.spoint))
end
struct PPExpand{N, P<:AbstractPID, D<:Number} <: VectorSpace{Operator{D, ID{OID{Index{P, NID{Char}}, SVector{N, D}}, 2}}}
    direction::SVector{N, D}
    points::NTuple{2, Point{N, P, D}}
end
@inline VectorSpaceStyle(::Type{<:PPExpand}) = VectorSpaceCartesian()
@inline shape(pnce::PPExpand) = (1:length(pnce.direction), 1:length(pnce.direction), 1:4)
function Operator(index::CartesianIndex{3}, pnce::PPExpand)
    dir₁ = Char(Int('x')+index[1]-1)
    dir₂ = Char(Int('x')+index[2]-1)
    coeff = index[3]∈(1, 4) ? 1 : -1
    pos₁, pos₂ = index[3]==1 ? (1, 1) : index[3]==2 ? (1, 2) : index[3]==3 ? (2, 1) : (2, 2)
    oid₁ = OID(Index(pnce.points[pos₁].pid, NID('u', dir₁)), pnce.points[pos₁].rcoord, pnce.points[pos₁].icoord)
    oid₂ = OID(Index(pnce.points[pos₂].pid, NID('u', dir₂)), pnce.points[pos₂].rcoord, pnce.points[pos₂].icoord)
    return Operator(pnce.direction[index[1]]*pnce.direction[index[2]]*coeff, ID(oid₁, oid₂))
end

"""
    kinetic"" -> Couplings

The kinetic energy part of phonon couplings.
"""
macro kinetic_str(::String) Couplings(PhononCoupling(1, ('p', 'p'))) end

"""
    potential"" -> Couplings

The potential energy part of phonon couplings.
"""
macro potential_str(::String) Couplings(PhononCoupling(1, ('u', 'u'))) end

"""
    PhononKinetic(id::Symbol, value::Any;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Kinetic energy part of phonons.

Type alias for `Term{:PhononKinetic, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const PhononKinetic{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PhononKinetic, id, V, Int, C, A, M}
@inline function PhononKinetic(id::Symbol, value::Any;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    Term{:PhononKinetic}(id, value, 0, couplings=kinetic"", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:PhononKinetic}) = :pnk
@inline ishermitian(::Type{<:PhononKinetic}) = true

"""
    PhononPotential(id::Symbol, value::Any, bondkind::Int;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Potential energy part of phonons.

Type alias for `Term{:PhononPotential, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`
"""
const PhononPotential{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PhononPotential, id, V, Int, C, A, M}
@inline function PhononPotential(id::Symbol, value::Any, bondkind::Int;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    Term{:PhononPotential}(id, value, bondkind, couplings=potential"", amplitude=amplitude, modulate=modulate)
end
@inline abbr(::Type{<:PhononPotential}) = :pnp
@inline ishermitian(::Type{<:PhononPotential}) = true
@inline couplingcenters(::PhononCoupling, ::Bond, ::Val{:PhononPotential}) = (1, 2)

"""
    PhononTerm

Type alias for `Union{PhononKinetic, PhononPotential}`.
"""
const PhononTerm = Union{PhononKinetic, PhononPotential}

end # module
