module SpinPackage

using StaticArrays: SVector
using Printf: @printf, @sprintf
using ..Spatials: PID, Point, Bond, AbstractBond
using ..DegreesOfFreedom: IID, Internal, Index, OID, AbstractCompositeOID, OIDToTuple, Operator, LaTeX, latexformat, Config, oidtype
using ..Terms: wildcard, Subscripts, SubID, subscriptsexpr, Coupling, Couplings, Term, TermCouplings, TermAmplitude, TermModulate
using ...Essentials: kind
using ...Prerequisites: Float, decimaltostr, delta
using ...Prerequisites.Traits: rawtype
using ...Mathematics.VectorSpaces: CartesianVectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID, ID

import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
import ..DegreesOfFreedom: script, latexname, isHermitian
import ..Terms: nonconstrain, couplingcenters, otype, abbr
import ...Interfaces: rank, expand, permute

export sdefaultlatex, usualspinindextotuple
export SID, Spin, SIndex, SOperator, SCID, SpinCoupling, SpinTerm, totalspin
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
    for i = 1:N, j = 1:N
        row, col = N+1-i, N+1-j
        m, n = S+1-i, S+1-j
        result[row, col] = (sid.tag == 'x') ? (delta(i+1, j)+delta(i, j+1))*sqrt(S*(S+1)-m*n)/2 :
            (sid.tag == 'y') ? (delta(i+1, j)-delta(i, j+1))*sqrt(S*(S+1)-m*n)/2im :
            (sid.tag == 'z') ? delta(i, j)*m :
            (sid.tag == '+') ? delta(i+1, j)*sqrt(S*(S+1)-m*n) :
            delta(i, j+1)*sqrt(S*(S+1)-m*n)
    end
    return result
end

"""
    Spin{S} <: Internal{SID{S}}

The spin interanl degrees of freedom.
"""
struct Spin{S} <: Internal{SID{S}}
    atom::Int
    norbital::Int
    function Spin{S}(atom::Int, norbital::Int) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
        new{S}(atom, norbital)
    end
end
@inline Base.Dims(sp::Spin) = (sp.norbital, length(sidtagmap))
@inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sid.orbital, sidseqmap[sid.tag])
@inline SID(index::CartesianIndex{2}, sp::Spin) = SID{totalspin(sp)}(index[1], sidtagmap[index[2]])
Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S

"""
    Spin{S}(; atom::Int=1, norbital::Int=1) where S

Construct a spin degrees of freedom.
"""
@inline Spin{S}(; atom::Int=1, norbital::Int=1) where S = Spin{S}(atom, norbital)

"""
    SIndex{S, P} <: Index{PID{P}, SID{S}}

The spin index.
"""
struct SIndex{S, P} <: Index{PID{P}, SID{S}}
    scope::P
    site::Int
    orbital::Int
    tag::Char
end
function SIndex{S}(scope::P, site::Int, orbital::Int, tag::Char) where {S, P}
    @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "SIndex error: not supported spin($S)."
    return SIndex{S, P}(scope, site, orbital, tag)
end
Base.show(io::IO, index::SIndex) = @printf io "SIndex{%s}(%s)" totalspin(index) join(repr.(values(index)), ", ")
@inline @generated function Base.replace(index::SIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){totalspin(index)}($(exprs...)))
end
@inline Base.union(::Type{P}, ::Type{I}) where {P<:PID, I<:SID} = SIndex{totalspin(I), fieldtype(P, :scope)}
@inline totalspin(index::SIndex) = totalspin(typeof(index))
@inline totalspin(::Type{<:SIndex{S}}) where S = S

"""
    SIndex(pid::PID, sid::SID) -> SIndex

Construct a spin index by a pid and an sid.
"""
@inline SIndex(pid::PID, sid::SID) = SIndex{totalspin(sid)}(values(pid)..., values(sid)...)

"""
    script(::Val{:site}, index::SIndex; kwargs...) -> Int
    script(::Val{:orbital}, index::SIndex; kwargs...) -> Int
    script(::Val{:tag}, index::SIndex; kwargs...) -> Char

Get the required script of a spin oid.
"""
@inline script(::Val{:site}, index::SIndex; kwargs...) = index.site
@inline script(::Val{:orbital}, index::SIndex; kwargs...) = index.orbital
@inline script(::Val{:tag}, index::SIndex; kwargs...) = index.tag

"""
    sdefaultlatex

The default LaTeX format for a spin oid.
"""
const soptdefaultlatex = LaTeX{(:tag,), (:site, :orbital)}('S')
@inline latexname(::Type{<:SIndex}) = Symbol("SIndex")
@inline latexname(::Type{<:OID{<:SIndex}}) = Symbol("OID{SIndex}")
latexformat(SIndex, soptdefaultlatex)
latexformat(OID{<:SIndex}, soptdefaultlatex)

"""
    usualspinindextotuple

Indicate that the choosed fields are `(:scope, :site, :orbital)` when converting a spin index to tuple.
"""
const usualspinindextotuple = OIDToTuple(:scope, :site, :orbital)

"""
    SOperator{V<:Number, I<:ID{AbstractCompositeOID{<:SIndex}}} <: Operator{V, I}

Spin operator.
"""
struct SOperator{V<:Number, I<:ID{AbstractCompositeOID{<:SIndex}}} <: Operator{V, I}
    value::V
    id::I
    SOperator(value::Number) = new{typeof(value), Tuple{}}(value, ())
    SOperator(value::Number, id::ID{AbstractCompositeOID{<:SIndex}}) = new{typeof(value), typeof(id)}(value, id)
end

"""
    permute(::Type{<:SOperator}, id₁::OID{<:SIndex}, id₂::OID{<:SIndex}) -> Tuple{Vararg{SOperator}}

Permute two fermionic oid and get the result.
"""
function permute(::Type{<:SOperator}, id₁::OID{<:SIndex}, id₂::OID{<:SIndex})
    @assert id₁.index≠id₂.index || id₁.rcoord≠id₂.rcoord || id₁.icoord≠id₂.icoord "permute error: permuted ids should not be equal to each other."
    if usualspinindextotuple(id₁.index)==usualspinindextotuple(id₂.index) && id₁.rcoord==id₂.rcoord && id₁.icoord==id₂.icoord
        @assert totalspin(id₁.index)==totalspin(id₂.index) "permute error: noncommutable ids should have the same spin field."
        if id₁.index.tag == 'x'
            id₂.index.tag=='y' && return (SOperator(+1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='z' && return (SOperator(-1im, ID(permutesoid(id₁, 'y'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='+' && return (SOperator(-1, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='-' && return (SOperator(+1, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
        elseif id₁.index.tag == 'y'
            id₂.index.tag=='x' && return (SOperator(-1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='z' && return (SOperator(+1im, ID(permutesoid(id₁, 'x'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='+' && return (SOperator(-1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='-' && return (SOperator(-1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
        elseif id₁.index.tag == 'z'
            id₂.index.tag=='x' && return (SOperator(+1im, ID(permutesoid(id₁, 'y'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='y' && return (SOperator(-1im, ID(permutesoid(id₁, 'x'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='+' && return (SOperator(+1, ID(id₂)), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='-' && return (SOperator(-1, ID(id₂)), SOperator(1, ID(id₂, id₁)))
        elseif id₁.index.tag == '+'
            id₂.index.tag=='x' && return (SOperator(+1, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='y' && return (SOperator(+1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='z' && return (SOperator(-1, ID(id₁)), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='-' && return (SOperator(+2, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
        elseif id₁.index.tag == '-'
            id₂.index.tag=='x' && return (SOperator(-1, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='y' && return (SOperator(1im, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='z' && return (SOperator(+1, ID(id₁)), SOperator(1, ID(id₂, id₁)))
            id₂.index.tag=='+' && return (SOperator(-2, ID(permutesoid(id₁, 'z'))), SOperator(1, ID(id₂, id₁)))
        end
    else
        return (SOperator(1, ID(id₂, id₁)),)
    end
end
@inline permutesoid(id::OID{<:SIndex}, tag::Char) = replace(id, index=replace(id.index, tag=tag))

"""
    SCID{T<:Tuple{Vararg{Char}}, A<:Tuple} <: SimpleID

The id of the atoms and tags part of a spin coupling.
"""
struct SCID{T<:Tuple{Vararg{Char}}, A<:Tuple} <: SimpleID
    tags::T
    atoms::A
    function SCID(tags::NTuple{N, Char}, atoms::NTuple{N}) where N
        @assert mapreduce(∈(('x', 'y', 'z', '+', '-')), &, tags) "SCID error: not supported tags($tags)."
        new{typeof(tags), typeof(atoms)}(tags, atoms)
    end
end

"""
    SpinCoupling{V, T<:Tuple, A<:Tuple, O<:Subscripts, I<:Tuple{SCID, SubID}} <: Coupling{V, I}

Spin coupling.
"""
struct SpinCoupling{V, T<:Tuple, A<:Tuple, O<:Subscripts, I<:Tuple{SCID, SubID}} <: Coupling{V, I}
    value::V
    tags::T
    atoms::A
    orbitals::O
    function SpinCoupling(value::Number, tags::Tuple{Vararg{Char}}, atoms::Tuple, orbitals::Subscripts)
        @assert length(tags)==length(atoms)==length(orbitals) "SpinCoupling error: dismatched tags, atoms and orbitals."
        scid, obid = SCID(tags, atoms), SubID(orbitals)
        new{typeof(value), typeof(tags), typeof(atoms), typeof(orbitals), Tuple{typeof(scid), typeof(obid)}}(value, tags, atoms, orbitals)
    end
end
@inline parameternames(::Type{<:SpinCoupling}) = (:value, :tags, :atoms, :orbitals, :id)
@inline isparameterbound(::Type{<:SpinCoupling}, ::Val{:tags}, ::Type{T}) where {T<:Tuple} = !isconcretetype(T)
@inline isparameterbound(::Type{<:SpinCoupling}, ::Val{:atoms}, ::Type{A}) where {A<:Tuple} = !isconcretetype(A)
@inline isparameterbound(::Type{<:SpinCoupling}, ::Val{:orbitals}, ::Type{O}) where {O<:Subscripts} = !isconcretetype(O)
@inline contentnames(::Type{<:SpinCoupling}) = (:value, :id, :orbitals)
@inline getcontent(sc::SpinCoupling, ::Val{:id}) = ID(SCID(sc.tags, sc.atoms), SubID(sc.orbitals))
@inline function SpinCoupling(value::Number, id::Tuple{SCID, SubID}, orbitals::Subscripts)
    SpinCoupling(value, id[1].tags, id[1].atoms, orbitals)
end
@inline rank(::Type{<:SpinCoupling{V, T} where V}) where T = fieldcount(T)

"""
    SpinCoupling(value::Number,
        tags::NTuple{N, Char};
        atoms::Union{NTuple{N, Int}, Nothing}=nothing,
        orbitals::Union{NTuple{N, Int}, Subscripts, Nothing}=nothing
        ) where N

Spin coupling.
"""
function SpinCoupling(value::Number,
        tags::NTuple{N, Char};
        atoms::Union{NTuple{N, Int}, Nothing}=nothing,
        orbitals::Union{Int, NTuple{N, Int}, Subscripts}=N
        ) where N
    isnothing(atoms) && (atoms = ntuple(i->wildcard, Val(N)))
    isa(orbitals, Subscripts) || (orbitals = Subscripts(orbitals))
    return SpinCoupling(value, tags, atoms, orbitals)
end

"""
    show(io::IO, sc::SpinCoupling)

Show a spin coupling.
"""
function Base.show(io::IO, sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    @printf io ", tags=%s" join(NTuple{rank(sc), String}("S"*sidrepmap[tag] for tag in sc.tags), "")
    any(sc.atoms .≠ wildcard) && @printf io ", atoms=[%s]" join(sc.atoms, " ")
    any(sc.orbitals .≠ wildcard) && @printf io ", orbitals=%s" string(sc.orbitals)
    @printf io ")"
end

"""
    repr(sc::SpinCoupling) -> String

Get the repr representation of a spin coupling.
"""
function Base.repr(sc::SpinCoupling)
    contents = String[]
    any(sc.atoms .≠ wildcard) && push!(contents, @sprintf "sl[%s]" join(sc.atoms, " "))
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
    return SpinCoupling(sc₁.value*sc₂.value, (sc₁.tags..., sc₂.tags...), (sc₁.atoms..., sc₁.atoms...), sc₁.orbitals*sc₂.orbitals)
end

"""
    expand(sc::SpinCoupling, points::NTuple{R, Point}, spins::NTuple{R, Spin}, ::Val) where R  -> Union{SCExpand, Tuple{}}

Expand a spin coupling with the given set of points and spin degrees of freedom.
"""
function expand(sc::SpinCoupling, points::NTuple{R, Point}, spins::NTuple{R, Spin}, ::Val) where R
    @assert rank(sc)==R "expand error: dismatched rank."
    for (i, atom) in enumerate(sc.atoms)
        isa(atom, Int) && atom≠spins[i].atom && return ()
    end
    obexpands = collect(expand(sc.orbitals, NTuple{rank(sc), Int}(spins[i].norbital for i = 1:rank(sc))))
    return SCExpand{totalspins(spins)}(sc.value, points, obexpands, sc.tags)
end
@generated totalspins(spins::NTuple{R, Spin}) where R = Tuple(totalspin(fieldtype(spins, i)) for i = 1:R)
struct SCExpand{SPS, V, N, D, P} <: CartesianVectorSpace{Tuple{V, ID{OID{SIndex, SVector{D, Float}}, N}}}
    value::V
    points::NTuple{N, Point{D, PID{P}}}
    obexpands::Vector{NTuple{N, Int}}
    tags::NTuple{N, Char}
    function SCExpand{SPS}(value::V, points::NTuple{N, Point{D, PID{P}}}, obexpands::Vector{NTuple{N, Int}}, tags::NTuple{N, Char}) where {SPS, V, N, D, P}
        return new{SPS, V, N, D, P}(value, points, obexpands, tags)
    end
end
@inline @generated function Base.eltype(::Type{SCExpand{SPS, V, N, D, P}}) where {SPS, V, N, D, P}
    return Tuple{V, Tuple{[OID{SIndex{SPS[i], P}, SVector{D, Float}} for i = 1:N]...}}
end
@inline Base.Dims(sce::SCExpand) = (length(sce.obexpands),)
@generated function Tuple(index::CartesianIndex{1}, sce::SCExpand{SPS, V, N}) where {SPS, V, N}
    exprs = []
    for i = 1:N
        spin = SPS[i]
        push!(exprs, quote
            pid, rcoord, icoord = sce.points[$i].pid, sce.points[$i].rcoord, sce.points[$i].icoord
            sid = SID{$spin}(sce.obexpands[index[1]][$i], sce.tags[$i])
            OID(SIndex(pid, sid), rcoord, icoord)
        end)
    end
    return Expr(:tuple, :(sce.value), Expr(:tuple, exprs...))
end

"""
    Heisenberg(mode::String="+-z"; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2) -> Couplings

The Heisenberg couplings.
"""
function Heisenberg(mode::String="+-z"; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2)
    @assert mode=="+-z" || mode=="xyz" "Heisenberg error: not supported mode($mode)."
    if mode == "+-z"
        sc₁ = SpinCoupling(1//2, ('+', '-'), atoms=atoms, orbitals=orbitals)
        sc₂ = SpinCoupling(1//2, ('-', '+'), atoms=atoms, orbitals=orbitals)
        sc₃ = SpinCoupling(1//1, ('z', 'z'), atoms=atoms, orbitals=orbitals)
    else
        sc₁ = SpinCoupling(1, ('x', 'x'), atoms=atoms, orbitals=orbitals)
        sc₂ = SpinCoupling(1, ('y', 'y'), atoms=atoms, orbitals=orbitals)
        sc₃ = SpinCoupling(1, ('z', 'z'), atoms=atoms, orbitals=orbitals)
    end
    return Couplings(sc₁, sc₂, sc₃)
end

"""
    Ising(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2) -> Couplings

The Ising couplings.
"""
@inline function Ising(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2)
    @assert tag in ('x', 'y', 'z') "Ising error: not supported input tag($tag)."
    return Couplings(SpinCoupling(1, (tag, tag), atoms=atoms, orbitals=orbitals))
end

"""
    Gamma(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2) -> Couplings

The Gamma couplings.
"""
function Gamma(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2)
    @assert tag in ('x', 'y', 'z') "Gamma error: not supported input tag($tag)."
    t₁, t₂ = tag=='x' ? ('y', 'z') : tag=='y' ? ('z', 'x') : ('x', 'y')
    sc₁ = SpinCoupling(1, (t₁, t₂), atoms=atoms, orbitals=orbitals)
    sc₂ = SpinCoupling(1, (t₂, t₁), atoms=atoms, orbitals=orbitals)
    return Couplings(sc₁, sc₂)
end

"""
    DM(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2) -> Couplings

The DM couplings.
"""
function DM(tag::Char; atoms::Union{NTuple{2, Int}, Nothing}=nothing, orbitals::Union{Int, NTuple{2, Int}, Subscripts}=2)
    @assert tag in ('x', 'y', 'z') "DM error: not supported input tag($tag)."
    t₁, t₂ = tag=='x' ? ('y', 'z') : tag=='y' ? ('z', 'x') : ('x', 'y')
    sc₁ = SpinCoupling(+1, (t₁, t₂), atoms=atoms, orbitals=orbitals)
    sc₂ = SpinCoupling(-1, (t₂, t₁), atoms=atoms, orbitals=orbitals)
    return Couplings(sc₁, sc₂)
end

@inline atomwrapper(::Nothing) = nothing
@inline atomwrapper(value::Int) = (value,)
@inline orbitalwrapper(::Nothing) = 1
@inline orbitalwrapper(value::Int) = (value,)
"""
    Sˣ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing) -> Couplings

The single Sˣ coupling.
"""
@inline function Sˣ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing)
    Couplings(SpinCoupling(1, ('x',), atoms=atomwrapper(atom), orbitals=orbitalwrapper(orbital)))
end

"""
    Sʸ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing) -> Couplings

The single Sʸ coupling.
"""
@inline function Sʸ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing)
    Couplings(SpinCoupling(1, ('y',), atoms=atomwrapper(atom), orbitals=orbitalwrapper(orbital)))
end

"""
    Sᶻ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing) -> Couplings

The single Sᶻ coupling.
"""
@inline function Sᶻ(; atom::Union{Int, Nothing}=nothing, orbital::Union{Int, Nothing}=nothing)
    Couplings(SpinCoupling(1, ('z',), atoms=atomwrapper(atom), orbitals=orbitalwrapper(orbital)))
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
    return SpinCoupling(coeff, tags; attrpairs...)
end
function sccomponent(str::AbstractString)
    attrname = str[1:2]=="sl" ? :atoms : str[1:2]=="ob" ? :orbitals : error("sccomponent error: wrong input pattern.")
    expr = Meta.parse(str[3:end])
    if attrname == :atoms
        @assert expr.head∈(:hcat, :vect) "sccomponent error: wrong input pattern for atoms."
        attrvalue = Tuple(expr.args)
    else
        @assert expr.head∈(:call, :hcat, :vcat, :vect) "sccomponent error: wrong input pattern for orbitals."
        attrvalue = eval(subscriptsexpr(expr))
    end
    return attrname => attrvalue
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
    heisenberg"sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    heisenberg"xyz sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    heisenberg"+-z sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings

The Heisenberg couplings.
"""
macro heisenberg_str(str::String)
    slice = findfirst(r"xyz|\+\-z", str)
    mode = isnothing(slice) ? "+-z" : str[slice]
    isnothing(slice) || (str = strip(str[nextind(str, slice.stop):end]))
    attrpairs = scpairs(str)
    return Heisenberg(mode; attrpairs...)
end

"""
    ising"x sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    ising"y sl[a₁ a₂] ⊗ ob[o₁ o₂])" -> Couplings
    ising"z sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings

The Ising couplings.
"""
macro ising_str(str::String)
    @assert str[1] ∈ ('x', 'y', 'z') "@ising_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@ising_str error: wrong input pattern."
    attrpairs = scpairs(str[3:end])
    return Ising(str[1]; attrpairs...)
end

"""
    gamma"x sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    gamma"y sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    gamma"z sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings

The Gamma couplings.
"""
macro gamma_str(str::String)
    @assert str[1] in ('x', 'y', 'z') "@gamma_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@gamma_str error: wrong input pattern."
    attrpairs = scpairs(str[3:end])
    return Gamma(str[1]; attrpairs...)
end

"""
    dm"x sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    dm"y sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings
    dm"z sl[a₁ a₂] ⊗ ob[o₁ o₂]" -> Couplings

The DM couplings.
"""
macro dm_str(str::String)
    @assert str[1] in ('x', 'y', 'z') "@dm_str error: wrong input pattern."
    length(str)>1 && @assert str[2]==' ' "@dm_str error: wrong input pattern."
    attrpairs = scpairs(str[3:end])
    return DM(str[1]; attrpairs...)
end

"""
    sˣ"sl[a]⊗ob[o]" -> Couplings
    sʸ"sl[a]⊗ob[o]" -> Couplings
    sᶻ"sl[a]⊗ob[o]" -> Couplings

The single Sˣ/Sʸ/Sᶻ coupling.
"""
macro sˣ_str(str::String) Couplings(SpinCoupling(1, ('x',); scpairs(str)...)) end
macro sʸ_str(str::String) Couplings(SpinCoupling(1, ('y',); scpairs(str)...)) end
macro sᶻ_str(str::String) Couplings(SpinCoupling(1, ('z',); scpairs(str)...)) end

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

"""
    otype(T::Type{<:Term}, C::Type{<:Config{<:Spin}}, B::Type{<:AbstractBond})

Get the operator type of a spin term.
"""
@inline otype(T::Type{<:Term}, C::Type{<:Config{<:Spin}}, B::Type{<:AbstractBond}) = SOperator{valtype(T), ID{oidtype(valtype(C), eltype(B), kind(T)|>Val), rank(T)}}

end # module
