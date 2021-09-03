module SpinPackage

using StaticArrays: SVector
using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, Point, Bond, AbstractBond
using ..DegreesOfFreedom: SimpleIID, SimpleInternal, Index, OID, AbstractCompositeOID, OIDToTuple, Operator, LaTeX, latexformat, Hilbert
using ..Terms: wildcard, Subscripts, SubID, subscriptsexpr, Coupling, Couplings, couplingpoints, couplinginternals, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, decimaltostr, delta
using ...Prerequisites.Traits: rawtype
using ...Mathematics.VectorSpaces: CartesianVectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID, ID

import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: nonconstrain, couplingcenters, abbr
import ...Interfaces: rank, expand, permute
import ...Mathematics.VectorSpaces: shape, ndimshape
import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent

export sdefaultlatex, usualspinindextotuple
export SID, Spin, SCID, SpinCoupling, SpinTerm, totalspin
export @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

const sidtagmap = Dict(1=>'x', 2=>'y', 3=>'z', 4=>'+', 5=>'-')
const sidseqmap = Dict(v=>k for (k, v) in sidtagmap)
const sidajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')
const sidrepmap = Dict('x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ', '+'=>'⁺', '-'=>'⁻', '0'=>'⁰')
const sidreprevmap = Dict(v=>k for (k, v) in sidrepmap)

"""
    SID{S} <: SimpleIID

The spin id.
"""
struct SID{S} <: SimpleIID
    orbital::Int
    tag::Char
    function SID{S}(orbital::Int, tag::Char) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "SID error: not supported spin($S)."
        @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
        new{S}(orbital, tag)
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
    SID{S}(tag::Char; orbital::Int=1) where S

Create a spin id.
"""
@inline SID{S}(tag::Char; orbital::Int=1) where S = SID{S}(orbital, tag)

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
    Spin{S} <: SimpleInternal{SID{S}}

The spin interanl degrees of freedom.
"""
struct Spin{S} <: SimpleInternal{SID{S}}
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
const soptdefaultlatex = LaTeX{(:tag,), (:site, :orbital)}('S')
@inline latexname(::Type{<:Index{<:AbstractPID, <:SID}}) = Symbol("Index{AbstractPID, SID}")
@inline latexname(::Type{<:AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}}) = Symbol("AbstractCompositeOID{Index{AbstractPID, SID}}")
latexformat(Index{<:AbstractPID, <:SID}, soptdefaultlatex)
latexformat(AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}, soptdefaultlatex)

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
    SCID{T<:Tuple{Vararg{Char}}} <: SimpleID

The id of the tags part of a spin coupling.
"""
struct SCID{T<:Tuple{Vararg{Char}}} <: SimpleID
    tags::T
    function SCID(tags::NTuple{N, Char}) where N
        @assert mapreduce(∈(('x', 'y', 'z', '+', '-')), &, tags) "SCID error: not supported tags($tags)."
        new{typeof(tags)}(tags)
    end
end

"""
    SpinCoupling{V, T<:Tuple, O<:Subscripts, I<:Tuple{SCID, SubID}} <: Coupling{V, I}

Spin coupling.
"""
struct SpinCoupling{V, T<:Tuple, O<:Subscripts, I<:Tuple{SCID, SubID}} <: Coupling{V, I}
    value::V
    tags::T
    orbitals::O
    function SpinCoupling(value::Number, tags::Tuple{Vararg{Char}}, orbitals::Subscripts)
        @assert length(tags)==length(orbitals) "SpinCoupling error: dismatched tags and orbitals."
        scid, obid = SCID(tags), SubID(orbitals)
        new{typeof(value), typeof(tags), typeof(orbitals), Tuple{typeof(scid), typeof(obid)}}(value, tags, orbitals)
    end
end
@inline parameternames(::Type{<:SpinCoupling}) = (:value, :tags, :orbitals, :id)
@inline isparameterbound(::Type{<:SpinCoupling}, ::Val{:tags}, ::Type{T}) where {T<:Tuple} = !isconcretetype(T)
@inline isparameterbound(::Type{<:SpinCoupling}, ::Val{:orbitals}, ::Type{O}) where {O<:Subscripts} = !isconcretetype(O)
@inline contentnames(::Type{<:SpinCoupling}) = (:value, :id, :orbitals)
@inline getcontent(sc::SpinCoupling, ::Val{:id}) = ID(SCID(sc.tags), SubID(sc.orbitals))
@inline function SpinCoupling(value::Number, id::Tuple{SCID, SubID}, orbitals::Subscripts)
    SpinCoupling(value, id[1].tags, orbitals)
end
@inline rank(::Type{<:SpinCoupling{V, T} where V}) where T = fieldcount(T)
function Base.show(io::IO, sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    @printf io ", tags=%s" join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.tags), "")
    any(sc.orbitals .≠ wildcard) && @printf io ", orbitals=%s" string(sc.orbitals)
    @printf io ")"
end

"""
    SpinCoupling(value::Number,
        tags::NTuple{N, Char};
        orbitals::Union{NTuple{N, Int}, Subscripts}=Subscripts(N)
        ) where N

Spin coupling.
"""
function SpinCoupling(value::Number,
        tags::NTuple{N, Char};
        orbitals::Union{NTuple{N, Int}, Subscripts}=Subscripts(N)
        ) where N
    isa(orbitals, Subscripts) || (orbitals = Subscripts(orbitals))
    return SpinCoupling(value, tags, orbitals)
end

"""
    repr(sc::SpinCoupling) -> String

Get the repr representation of a spin coupling.
"""
function Base.repr(sc::SpinCoupling)
    contents = String[]
    any(sc.orbitals .≠ wildcard) && push!(contents, @sprintf "ob%s" sc.orbitals)
    result = @sprintf "%s %s" decimaltostr(sc.value) join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.tags), "")
    length(contents)>0 && (result = @sprintf "%s %s" result join(contents, " ⊗ "))
    return result
end

"""
    *(sc₁::SpinCoupling, sc₂::SpinCoupling) -> SpinCoupling

Get the multiplication between two spin couplings.
"""
@inline function Base.:*(sc₁::SpinCoupling, sc₂::SpinCoupling)
    return SpinCoupling(sc₁.value*sc₂.value, (sc₁.tags..., sc₂.tags...), sc₁.orbitals*sc₂.orbitals)
end

"""
    expand(sc::SpinCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val) -> SCExpand

Expand a spin coupling with the given set of points and Hilbert space.
"""
function expand(sc::SpinCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)
    points = couplingpoints(sc, bond, info)
    spins = couplinginternals(sc, bond, hilbert, info)
    @assert rank(sc)==length(points)==length(spins) "expand error: dismatched rank."
    obexpands = collect(expand(sc.orbitals, NTuple{rank(sc), Int}(spins[i].norbital for i = 1:rank(sc))))
    return SCExpand{totalspin(spins)}(sc.value, points, obexpands, sc.tags)
end
@generated totalspin(spins::NTuple{R, Spin}) where R = Tuple(totalspin(fieldtype(spins, i)) for i = 1:R)
struct SCExpand{SPS, V, N, D, P<:AbstractPID, DT<:Number} <: CartesianVectorSpace{Tuple{V, ID{OID{<:Index, SVector{D, DT}}, N}}}
    value::V
    points::NTuple{N, Point{D, P, DT}}
    obexpands::Vector{NTuple{N, Int}}
    tags::NTuple{N, Char}
    function SCExpand{SPS}(value::V, points::NTuple{N, Point{D, P, DT}}, obexpands::Vector{NTuple{N, Int}}, tags::NTuple{N, Char}) where {SPS, V, N, D, P<:AbstractPID, DT<:Number}
        return new{SPS, V, N, D, P, DT}(value, points, obexpands, tags)
    end
end
@inline @generated function Base.eltype(::Type{SCExpand{SPS, V, N, D, P, DT}}) where {SPS, V, N, D, P<:AbstractPID, DT<:Number}
    return Tuple{V, Tuple{[OID{Index{P, SID{SPS[i]}}, SVector{D, DT}} for i = 1:N]...}}
end
@inline shape(sce::SCExpand) = (1:length(sce.obexpands),)
@inline ndimshape(::Type{<:SCExpand}) = 1
@generated function Tuple(index::CartesianIndex{1}, sce::SCExpand{SPS, V, N}) where {SPS, V, N}
    exprs = []
    for i = 1:N
        spin = SPS[i]
        push!(exprs, quote
            pid, rcoord, icoord = sce.points[$i].pid, sce.points[$i].rcoord, sce.points[$i].icoord
            sid = SID{$spin}(sce.obexpands[index[1]][$i], sce.tags[$i])
            OID(Index(pid, sid), rcoord, icoord)
        end)
    end
    return Expr(:tuple, :(sce.value), Expr(:tuple, exprs...))
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
    length(str)==0 && return Expr(:kw, :orbitals, Subscripts(n))
    @assert str[1:2]=="ob" "scorbitals error: wrong input pattern."
    expr = Meta.parse(str[3:end])
    @assert expr.head∈(:call, :hcat, :vcat, :vect) "scorbitals error: wrong input pattern for orbitals."
    attrvalue = subscriptsexpr(expr)
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
function heisenbergxyz(;orbitals=Subscripts(2))
    return Couplings(
        SpinCoupling(1, ('x', 'x'), orbitals=orbitals),
        SpinCoupling(1, ('y', 'y'), orbitals=orbitals),
        SpinCoupling(1, ('z', 'z'), orbitals=orbitals)
    )
end
function heisenbergpmz(;orbitals=Subscripts(2))
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
function gamma(t₁::Char, t₂::Char; orbitals=Subscripts(2))
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
function dm(t₁::Char, t₂::Char; orbitals=Subscripts(2))
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

Type alias for `Term{:SpinTerm, R, id, V, <:Any, <:TermCouplings, <:TermAmplitude, <:TermModulate}`.
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
