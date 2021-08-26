module SpinPackage

using StaticArrays: SVector
using Printf: @printf, @sprintf
using ..Spatials: AbstractPID, Point, Bond, AbstractBond
using ..DegreesOfFreedom: IID, Internal, Index, OID, AbstractCompositeOID, OIDToTuple, Operator, LaTeX, latexformat, Config
using ..Terms: wildcard, Subscripts, SubID, subscriptsexpr, Coupling, Couplings, couplingpoints, couplinginternals, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, decimaltostr, delta
using ...Prerequisites.Traits: rawtype
using ...Mathematics.VectorSpaces: CartesianVectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID, ID

import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: nonconstrain, couplingcenters, abbr
import ...Interfaces: rank, expand, permute

export sdefaultlatex, usualspinindextotuple
export SID, Spin, SCID, SpinCoupling, SpinTerm, totalspin
export Heisenberg, Ising, Gamma, DM, Sˣ, Sʸ, Sᶻ
export @sc_str, @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str

const sidtagmap = Dict(1=>'x', 2=>'y', 3=>'z', 4=>'+', 5=>'-')
const sidseqmap = Dict(v=>k for (k, v) in sidtagmap)
const sidajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')
const sidrepmap = Dict('x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ', '+'=>'⁺', '-'=>'⁻', '0'=>'⁰')
const sidreprevmap = Dict(v=>k for (k, v) in sidrepmap)

"""
    SID{S} <: IID

The spin id.
"""
struct SID{S} <: IID
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
    Spin{S} <: Internal{SID{S}}

The spin interanl degrees of freedom.
"""
struct Spin{S} <: Internal{SID{S}}
    norbital::Int
    function Spin{S}(norbital::Int) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
        new{S}(norbital)
    end
end
@inline Base.Dims(sp::Spin) = (sp.norbital, length(sidtagmap))
@inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sid.orbital, sidseqmap[sid.tag])
@inline SID(index::CartesianIndex{2}, sp::Spin) = SID{totalspin(sp)}(index[1], sidtagmap[index[2]])
Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S

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

Indicate that the choosed fields are `(:scope, :site, :orbital)` when converting a spin index to tuple.
"""
const usualspinindextotuple = OIDToTuple(:scope, :site, :orbital)

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
    show(io::IO, sc::SpinCoupling)

Show a spin coupling.
"""
function Base.show(io::IO, sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    @printf io ", tags=%s" join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.tags), "")
    any(sc.orbitals .≠ wildcard) && @printf io ", orbitals=%s" string(sc.orbitals)
    @printf io ")"
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
    expand(sc::SpinCoupling, bond::AbstractBond, config::Config, info::Val) -> SCExpand

Expand a spin coupling with the given set of points and spin degrees of freedom.
"""
function expand(sc::SpinCoupling, bond::AbstractBond, config::Config, info::Val)
    points = couplingpoints(sc, bond, info)
    spins = couplinginternals(sc, bond, config, info)
    @assert rank(sc)==length(points)==length(spins) "expand error: dismatched rank."
    obexpands = collect(expand(sc.orbitals, NTuple{rank(sc), Int}(spins[i].norbital for i = 1:rank(sc))))
    return SCExpand{totalspins(spins)}(sc.value, points, obexpands, sc.tags)
end
@generated totalspins(spins::NTuple{R, Spin}) where R = Tuple(totalspin(fieldtype(spins, i)) for i = 1:R)
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
@inline Base.Dims(sce::SCExpand) = (length(sce.obexpands),)
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
    Heisenberg(mode::String="+-z"; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2)) -> Couplings

The Heisenberg couplings.
"""
function Heisenberg(mode::String="+-z"; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2))
    @assert mode=="+-z" || mode=="xyz" "Heisenberg error: not supported mode($mode)."
    if mode == "+-z"
        sc₁ = SpinCoupling(1//2, ('+', '-'), orbitals=orbitals)
        sc₂ = SpinCoupling(1//2, ('-', '+'), orbitals=orbitals)
        sc₃ = SpinCoupling(1//1, ('z', 'z'), orbitals=orbitals)
    else
        sc₁ = SpinCoupling(1, ('x', 'x'), orbitals=orbitals)
        sc₂ = SpinCoupling(1, ('y', 'y'), orbitals=orbitals)
        sc₃ = SpinCoupling(1, ('z', 'z'), orbitals=orbitals)
    end
    return Couplings(sc₁, sc₂, sc₃)
end

"""
    Ising(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2)) -> Couplings

The Ising couplings.
"""
@inline function Ising(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2))
    @assert tag in ('x', 'y', 'z') "Ising error: not supported input tag($tag)."
    return Couplings(SpinCoupling(1, (tag, tag), orbitals=orbitals))
end

"""
    Gamma(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2)) -> Couplings

The Gamma couplings.
"""
function Gamma(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2))
    @assert tag in ('x', 'y', 'z') "Gamma error: not supported input tag($tag)."
    t₁, t₂ = tag=='x' ? ('y', 'z') : tag=='y' ? ('z', 'x') : ('x', 'y')
    sc₁ = SpinCoupling(1, (t₁, t₂), orbitals=orbitals)
    sc₂ = SpinCoupling(1, (t₂, t₁), orbitals=orbitals)
    return Couplings(sc₁, sc₂)
end

"""
    DM(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2)) -> Couplings

The DM couplings.
"""
function DM(tag::Char; orbitals::Union{NTuple{2, Int}, Subscripts}=Subscripts(2))
    @assert tag in ('x', 'y', 'z') "DM error: not supported input tag($tag)."
    t₁, t₂ = tag=='x' ? ('y', 'z') : tag=='y' ? ('z', 'x') : ('x', 'y')
    sc₁ = SpinCoupling(+1, (t₁, t₂), orbitals=orbitals)
    sc₂ = SpinCoupling(-1, (t₂, t₁), orbitals=orbitals)
    return Couplings(sc₁, sc₂)
end

@inline orbitalwrapper(value::Int) = (value,)
@inline orbitalwrapper(::Symbol) = Subscripts(1)
"""
    Sˣ(; orbital::Union{Int, Symbol}=wildcard) -> Couplings

The single Sˣ coupling.
"""
@inline function Sˣ(; orbital::Union{Int, Symbol}=wildcard)
    Couplings(SpinCoupling(1, ('x',), orbitals=orbitalwrapper(orbital)))
end

"""
    Sʸ(; orbital::Union{Int, Symbol}=wildcard) -> Couplings

The single Sʸ coupling.
"""
@inline function Sʸ(; orbital::Union{Int, Symbol}=wildcard)
    Couplings(SpinCoupling(1, ('y',), orbitals=orbitalwrapper(orbital)))
end

"""
    Sᶻ(; orbital::Union{Int, Symbol}=wildcard) -> Couplings

The single Sᶻ coupling.
"""
@inline function Sᶻ(; orbital::Union{Int, Symbol}=wildcard)
    Couplings(SpinCoupling(1, ('z',), orbitals=orbitalwrapper(orbital)))
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
    attrpairs = scpairs(strip(str[nextind(str, lpos):end]))
    return Expr(:call, :SpinCoupling, Expr(:parameters, attrpairs...), coeff, tags)
end
function sccomponent(str::AbstractString)
    @assert str[1:2]=="ob" "sccomponent error: wrong input pattern."
    expr = Meta.parse(str[3:end])
    @assert expr.head∈(:call, :hcat, :vcat, :vect) "sccomponent error: wrong input pattern for orbitals."
    attrvalue = subscriptsexpr(expr)
    return Expr(:kw, :orbitals, attrvalue)
end
function scpairs(str::AbstractString)
    attrpairs = []
    if length(str) > 0
        for component in split(replace(str, " ⊗ "=>"⊗"), '⊗')
            push!(attrpairs, sccomponent(component))
        end
    end
    return attrpairs
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
    attrpairs = scpairs(str)
    return Expr(:call, :Heisenberg, Expr(:parameters, attrpairs...), mode)
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
    attrpairs = scpairs(str[3:end])
    return Expr(:call, :Ising, Expr(:parameters, attrpairs...), str[1])
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
    attrpairs = scpairs(str[3:end])
    return Expr(:call, :Gamma, Expr(:parameters, attrpairs...), str[1])
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
    attrpairs = scpairs(str[3:end])
    return Expr(:call, :DM, Expr(:parameters, attrpairs...), str[1])
end

"""
    sˣ"ob[o]" -> Couplings
    sʸ"ob[o]" -> Couplings
    sᶻ"ob[o]" -> Couplings

The single Sˣ/Sʸ/Sᶻ coupling.
"""
macro sˣ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, Expr(:parameters, scpairs(str)...), 1, ('x',))) end
macro sʸ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, Expr(:parameters, scpairs(str)...), 1, ('y',))) end
macro sᶻ_str(str::String) Expr(:call, :Couplings, Expr(:call, :SpinCoupling, Expr(:parameters, scpairs(str)...), 1, ('z',))) end

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
