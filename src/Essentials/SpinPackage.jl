module SpinPackage

using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, Bond
using ..DegreesOfFreedom: SimpleIID, SimpleInternal, Index, OID, AbstractCompositeOID, OIDToTuple, Operator, LaTeX, latexformat
using ..Terms: IIDSpace, IIDConstrain, ConstrainID, wildcard, Subscript, subscriptexpr, diagonal
using ..Terms: Coupling, Couplings, Term, TermCouplings, TermAmplitude, TermModulate
using ...Interfaces: rank
using ...Prerequisites: Float, decimaltostr, delta
using ...Prerequisites.Traits: rawtype
using ...Mathematics.AlgebraOverFields: ID

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: couplingcenters, abbr
import ...Interfaces: permute
import ...Mathematics.VectorSpaces: shape, ndimshape

export sdefaultlatex, usualspinindextotuple
export SID, Spin, SpinCoupling, SpinTerm, totalspin
export @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

const sidtagmap = Dict(1=>'x', 2=>'y', 3=>'z', 4=>'+', 5=>'-')
const sidseqmap = Dict(v=>k for (k, v) in sidtagmap)
const sidajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')
const sidrepmap = Dict('x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ', '+'=>'⁺', '-'=>'⁻', '0'=>'⁰')
const sidreprevmap = Dict(v=>k for (k, v) in sidrepmap)

"""
    SID{S, O<:Union{Int, Symbol}} <: SimpleIID

The spin id.
"""
struct SID{S, O<:Union{Int, Symbol}} <: SimpleIID
    orbital::O
    tag::Char
    function SID{S}(orbital::Union{Int, Symbol}, tag::Char) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) || S==wildcard "SID error: not supported spin($S)."
        @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
        new{S, typeof(orbital)}(orbital, tag)
    end
end
@inline Base.adjoint(sid::SID) = SID{totalspin(sid)}(sid.orbital, sidajointmap[sid.tag])
Base.show(io::IO, sid::SID) = @printf io "SID{%s}(%s)" totalspin(sid) join(repr.(values(sid)), ", ")
@inline @generated function Base.replace(sid::SID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(sid, $name))) for name in QuoteNode.(fieldnames(sid))]
    return :(rawtype(typeof(sid)){totalspin(sid)}($(exprs...)))
end
@inline totalspin(sid::SID) = totalspin(typeof(sid))
@inline totalspin(::Type{<:SID{S}}) where S = S

"""
    SID{S}(tag::Char; orbital::Union{Int, Symbol}=1) where S

Create a spin id.
"""
@inline SID{S}(tag::Char; orbital::Union{Int, Symbol}=1) where S = SID{S}(orbital, tag)

"""
    Matrix(sid::SID{S}, dtype::Type{<:Number}=Complex{Float}) where S -> Matrix{dtype}

Get the matrix representation of a sid.
"""
function Base.Matrix(sid::SID{S}, dtype::Type{<:Number}=Complex{Float}) where S
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
    Spin{S} <: SimpleInternal{SID{S, Int}}

The spin interanl degrees of freedom.
"""
struct Spin{S} <: SimpleInternal{SID{S, Int}}
    norbital::Int
    function Spin{S}(norbital::Int) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
        new{S}(norbital)
    end
end
@inline shape(sp::Spin) = (1:sp.norbital, 1:length(sidtagmap))
@inline ndimshape(::Type{<:Spin}) = 2
@inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sid.orbital, sidseqmap[sid.tag])
@inline SID(index::CartesianIndex{2}, sp::Spin) = SID{totalspin(sp)}(index[1], sidtagmap[index[2]])
Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S
function Base.show(io::IO, spin::Spin)
    @printf io "%s{%s}(%s)" spin|>typeof|>nameof totalspin(spin) join(("$name=$(getfield(spin, name))" for name in spin|>typeof|>fieldnames), ", ")
end
@inline function shape(iidspace::IIDSpace{<:Union{SID{wildcard, Symbol}, SID{S, Symbol}}, Spin{S}}) where S
    norbital, tag = iidspace.internal.norbital, sidseqmap[iidspace.iid.tag]
    return (1:norbital, tag:tag)
end
@inline function shape(iidspace::IIDSpace{<:Union{SID{wildcard, Int}, SID{S, Int}}, Spin{S}}) where S
    orbital, tag = iidspace.iid.orbital, sidseqmap[iidspace.iid.tag]
    @assert orbital<iidspace.internal.norbital+1 "shape error: orbital out of range."
    return (orbital:orbital, tag:tag)
end

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
    sdefaultlatex

The default LaTeX format for a spin oid.
"""
const sdefaultlatex = LaTeX{(:tag,), (:site,)}('S')
@inline latexname(::Type{<:Index{<:AbstractPID, <:SID}}) = Symbol("Index{AbstractPID, SID}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, SID}}")
latexformat(Index{<:AbstractPID, <:SID}, sdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}, sdefaultlatex)

"""
    usualspinindextotuple

Indicate that the choosed fields are `(:site, :orbital)` when converting a spin index to tuple.
"""
const usualspinindextotuple = OIDToTuple(:site, :orbital)

"""
    permute(id₁::OID{<:Index{<:AbstractPID, <:SID}}, id₂::OID{<:Index{<:AbstractPID, <:SID}}) -> Tuple{Vararg{Operator}}

Permute two spin oid and get the result.
"""
function permute(id₁::OID{<:Index{<:AbstractPID, <:SID}}, id₂::OID{<:Index{<:AbstractPID, <:SID}})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if usualspinindextotuple(id₁.index)==usualspinindextotuple(id₂.index) && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        @assert totalspin(id₁.index.iid)==totalspin(id₂.index.iid) "permute error: noncommutable ids should have the same spin field."
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

Type alias for "Coupling{V, I<:ID{SID}, C<:IIDConstrain, CI<:ConstrainID}".
"""
const SpinCoupling{V, I<:ID{SID}, C<:IIDConstrain, CI<:ConstrainID} = Coupling{V, I, C, CI}
@inline function SpinCoupling(value::Number, tags::NTuple{N, Char}, orbitals::Subscript{<:NTuple{N, Union{Int, Symbol}}}) where N
    return Coupling(value, ID(SID{wildcard}, orbitals.pattern, tags), IIDConstrain((orbital=orbitals,)))
end
function Base.show(io::IO, sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    @printf io ", tags=%s" join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.cid.tags), "")
    sc.cid.orbitals≠ntuple(i->wildcard, Val(rank(sc))) && @printf io ", orbitals=%s" repr(sc.constrain, 1:length(sc.constrain), :orbital)
    @printf io ")"
end
function Base.repr(sc::SpinCoupling)
    result = [@sprintf "%s %s" decimaltostr(sc.value) join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.cid.tags), "")]
    sc.cid.orbitals≠ntuple(i->wildcard, Val(rank(sc))) && push!(result, @sprintf "ob%s" repr(sc.constrain, 1:length(sc.constrain), :orbital))
    return join(result, " ")
end

"""
    SpinCoupling(value::Number, tags::NTuple{N, Char}; orbitals::Union{NTuple{N, Int}, Subscript}=Subscript(N)) where N

Spin coupling.
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
    SpinTerm{R}(id::Symbol, value::Any, bondkind::Any;
        couplings::Union{Function, Coupling, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false,
        ) where R

Spin term.

Type alias for `Term{:SpinTerm, R, id, V, B<:Any, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinTerm{R, id, V, B<:Any, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate} = Term{:SpinTerm, R, id, V, B, C, A, M}
@inline function SpinTerm{R}(id::Symbol, value::Any, bondkind::Any;
        couplings::Union{Function, Coupling, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        ) where R
    isa(couplings, TermCouplings) || (couplings = TermCouplings(couplings))
    isa(amplitude, TermAmplitude) || (amplitude = TermAmplitude(amplitude))
    isa(modulate, TermModulate) || (modulate = TermModulate(id, modulate))
    Term{:SpinTerm, R, id}(value, bondkind, couplings, amplitude, modulate, 1)
end
@inline abbr(::Type{<:SpinTerm}) = :sp
@inline isHermitian(::Type{<:SpinTerm}) = true
@inline function couplingcenters(sc::SpinCoupling, ::Bond, ::Val{:SpinTerm})
    @assert rank(sc)%2==0 "couplingcenters error: the rank of the input spin coupling should be even."
    return ntuple(i->2-i%2, Val(rank(sc)))
end

end # module
