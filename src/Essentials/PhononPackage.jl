module PhononPackage

using LinearAlgebra: norm
using StaticArrays: SVector
using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, Point, Bond, rcoord
using ..DegreesOfFreedom: SimpleIID, SimpleInternal, Index, AbstractCompositeOID, LaTeX, latexformat, OIDToTuple, OID, Operator, Hilbert
using ..Terms: IIDSpace, IIDConstrain, ConstrainID, Subscript, wildcard
using ..Terms: Coupling, Couplings, couplinginternals, Term, TermCouplings, TermAmplitude, TermModulate
using ...Interfaces: rank
using ...Prerequisites: atol, rtol, decimaltostr
using ...Prerequisites.Traits: getcontent
using ...Mathematics.AlgebraOverFields: ID
using ...Mathematics.VectorSpaces: CartesianVectorSpace

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: couplingcenters, abbr
import ...Interfaces: permute, expand
import ...Mathematics.VectorSpaces: shape, ndimshape

export pndefaultlatex, usualphononindextotuple
export PNID, Phonon, PhononCoupling, PhononKinetic, PhononPotential
export @kinetic_str, @potential_str

"""
    PNID{D<:Union{Char, Symbol}} <: SimpleIID

The phonon id.
"""
struct PNID{D<:Union{Char, Symbol}} <: SimpleIID
    tag::Char
    dir::D
    function PNID(tag::Char, dir::Union{Char, Symbol}=wildcard)
        @assert tag∈('p', 'u') "PNID error: wrong tag($tag)."
        isa(dir, Char) && @assert dir∈('x', 'y', 'z') "PNID error: wrong direction($dir)."
        new{typeof(dir)}(tag, dir)
    end
end
@inline Base.adjoint(pnid::PNID) = pnid

"""
    Phonon <: SimpleInternal{PNID{Char}}

The phonon internal degrees of freedom.
"""
struct Phonon <: SimpleInternal{PNID{Char}}
    ndir::Int
    function Phonon(ndir::Integer)
        @assert ndir∈(1, 2, 3) "Phonon error: wrong number of directions."
        new(ndir)
    end
end
@inline shape(pn::Phonon) = (1:2, 1:pn.ndir)
@inline ndimshape(::Type{Phonon}) = 2
@inline Base.CartesianIndex(pnid::PNID{Char}, ::Phonon) = CartesianIndex(pnid.tag=='u' ? 1 : 2, Int(pnid.dir)-Int('x')+1)
@inline PNID(index::CartesianIndex{2}, ::Phonon) = PNID(index[1]==1 ? 'u' : 'p', Char(Int('x')+index[2]-1))
@inline shape(iidspace::IIDSpace{PNID{Symbol}, Phonon}) = (iidspace.iid.tag=='u' ? (1:1) : (2:2), 1:iidspace.internal.ndir)
@inline function shape(iidspace::IIDSpace{PNID{Char}, Phonon})
    dir = Int(iidspace.iid.dir)-Int('x')+1
    @assert 0<dir<iidspace.internal.ndir+1 "shape error: dir out of range."
    return (iidspace.iid.tag=='u' ? (1:1) : (2:2), dir:dir)
end

"""
    script(::Val{:BD}, index::Index{<:AbstractPID, <:PNID}, l::LaTeX) -> Char
    script(::Val{:BD}, oid::AbstractCompositeOID{<:Index{<:AbstractPID, <:PNID}}, l::LaTeX) -> Char
    script(::Val{:site}, index::Index{<:AbstractPID, <:PNID}; kwargs...) -> Int
    script(::Val{:dir}, index::Index{<:AbstractPID, <:PNID}; kwargs...) -> Char

Get the required script of a phonon oid.
"""
@inline script(::Val{:BD}, index::Index{<:AbstractPID, <:PNID}, l::LaTeX) = l.body[index.iid.tag]
@inline script(::Val{:BD}, oid::AbstractCompositeOID{<:Index{<:AbstractPID, <:PNID}}, l::LaTeX) = l.body[getcontent(oid, :index).iid.tag]
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:PNID}; kwargs...) = index.pid.site
@inline script(::Val{:dir}, index::Index{<:AbstractPID, <:PNID}; kwargs...) = index.iid.dir

"""
    pndefaultlatex

The default LaTeX format for a phonon oid.
"""
const pndefaultlatex = LaTeX{(), (:site, :dir)}(Dict('p'=>'p', 'u'=>'u'), "", "")
@inline latexname(::Type{<:Index{<:AbstractPID, <:PNID}}) = Symbol("Index{AbstractPID, PNID}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, <:PNID}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, PNID}}")
latexformat(Index{<:AbstractPID, <:PNID}, pndefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, <:PNID}}, pndefaultlatex)

"""
    usualphononindextotuple

Indicate that the choosed fields are `(:tag, :site, :dir)` when converting a phonon index to tuple.
"""
const usualphononindextotuple = OIDToTuple(:tag, :site, :dir)

"""
    permute(id₁::OID{<:Index{<:AbstractPID, PNID{Char}}}, id₂::OID{<:Index{<:AbstractPID, PNID{Char}}}) -> Tuple{Vararg{Operator}}

Permute two phonon oids and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, PNID{Char}}}, id₂::OID{<:Index{<:AbstractPID, PNID{Char}}})
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

Type alias for `Coupling{V<:Number, I<:ID{PNID}, C<:IIDConstrain, CI<:ConstrainID}`.
"""
const PhononCoupling{V<:Number, I<:ID{PNID}, C<:IIDConstrain, CI<:ConstrainID} = Coupling{V, I, C, CI}
@inline function PhononCoupling(value::Number, tags::NTuple{N, Char}, dirs::Subscript{<:Union{NTuple{N, Char}, NTuple{N, Symbol}}}) where N
    return Coupling(value, ID(PNID, tags, dirs.pattern), IIDConstrain((dir=dirs,)))
end
function Base.show(io::IO, pnc::PhononCoupling)
    @printf io "PhononCoupling(value=%s" decimaltostr(pnc.value)
    @printf io ", tags=[%s]" join(pnc.cid.tags, " ")
    pnc.cid.dirs≠ntuple(i->wildcard, Val(rank(pnc))) && @printf io ", dirs=%s" repr(pnc.constrain, 1:length(pnc.constrain), :dir)
    @printf io ")"
end
function Base.repr(pnc::PhononCoupling)
    result = [@sprintf "%s [%s]" decimaltostr(pnc.value) join(pnc.cid.tags, " ")]
    pnc.cid.dirs≠ntuple(i->wildcard, Val(rank(pnc))) && push!(result, @sprintf "dr%s" repr(pnc.constrain, 1:length(pnc.constrain), :dir))
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
    expand(pnc::PhononCoupling{<:Number, <:ID{PNID{Symbol}}}, bond::Bond, hilbert::Hilbert, info::Val{:PhononPotential}) -> PPExpand

Expand the default phonon potential coupling on a given bond.
"""
function expand(pnc::PhononCoupling{<:Number, <:ID{PNID{Symbol}}}, bond::Bond, hilbert::Hilbert, info::Val{:PhononPotential})
    R̂ = rcoord(bond)/norm(rcoord(bond))
    @assert pnc.cid.tags==('u', 'u') "expand error: wrong tags of phonon coupling."
    @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of phonon coupling."
    pn₁, pn₂ = couplinginternals(pnc, bond, hilbert, info)
    @assert pn₁.ndir==pn₂.ndir==length(R̂) "expand error: dismatched number of directions."
    return PPExpand(R̂, (bond.epoint, bond.spoint))
end
struct PPExpand{N, P<:AbstractPID, D<:Number} <: CartesianVectorSpace{Tuple{D, ID{OID{Index{P, PNID{Char}}, SVector{N, D}}, 2}}}
    direction::SVector{N, D}
    points::NTuple{2, Point{N, P, D}}
end
@inline shape(pnce::PPExpand) = (1:length(pnce.direction), 1:3)
function Tuple(index::CartesianIndex{2}, pnce::PPExpand)
    dir = Char(Int('x')+index[1]-1)
    coeff = index[2]==2 ? -2 : 1
    pos₁, pos₂ = index[2]==1 ? (1, 1) : index[2]==2 ? (1, 2) : (2, 2)
    oid₁ = OID(Index(pnce.points[pos₁].pid, PNID('u', dir)), pnce.points[pos₁].rcoord, pnce.points[pos₁].icoord)
    oid₂ = OID(Index(pnce.points[pos₂].pid, PNID('u', dir)), pnce.points[pos₂].rcoord, pnce.points[pos₂].icoord)
    return (pnce.direction[index[1]])^2*coeff, ID(oid₁, oid₂)
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

Type alias for `Term{:PhononKinetic, 2, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const PhononKinetic{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PhononKinetic, 2, id, V, Int, C, A, M}
@inline function PhononKinetic(id::Symbol, value::Any;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    isa(amplitude, TermAmplitude) || (amplitude = TermAmplitude(amplitude))
    isa(modulate, TermModulate) || (modulate = TermModulate(id, modulate))
    Term{:PhononKinetic, 2, id}(value, 0, TermCouplings(kinetic""), amplitude, modulate, 1)
end
@inline abbr(::Type{<:PhononKinetic}) = :pnk
@inline isHermitian(::Type{<:PhononKinetic}) = true

"""
    PhononPotential(id::Symbol, value::Any, bondkind::Int;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )

Potential energy part of phonons.

Type alias for `Term{:PhononPotential, 2, id, V, Int, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`
"""
const PhononPotential{id, V, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:PhononPotential, 2, id, V, Int, C, A, M}
@inline function PhononPotential(id::Symbol, value::Any, bondkind::Int;
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        )
    isa(amplitude, TermAmplitude) || (amplitude = TermAmplitude(amplitude))
    isa(modulate, TermModulate) || (modulate = TermModulate(id, modulate))
    Term{:PhononPotential, 2, id}(value, bondkind, TermCouplings(potential""), amplitude, modulate, 1)
end
@inline abbr(::Type{<:PhononPotential}) = :pnp
@inline isHermitian(::Type{<:PhononPotential}) = true
@inline couplingcenters(::PhononCoupling, ::Bond, ::Val{:PhononPotential}) = (1, 2)

end # module
