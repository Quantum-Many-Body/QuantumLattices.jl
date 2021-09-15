using Test
using StaticArrays: SVector
using LinearAlgebra: dot
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.QuantumAlgebras: ID, sequence
using QuantumLattices.Essentials.Spatials: AbstractPID, PID, CPID, Point, pidtype, rcoord, icoord
using QuantumLattices.Essentials: update!, reset!
using QuantumLattices.Interfaces: decompose, rank, ⊗
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Traits: contentnames

import QuantumLattices.Essentials.DegreesOfFreedom: statistics, latexname, script
import QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

struct DID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.adjoint(sl::DID{Int}) = DID(3-sl.nambu)
@inline statistics(::Type{<:DID}) = :f

struct DFock <: SimpleInternal{DID{Int}}
    nnambu::Int
end
@inline shape(f::DFock) = (1:f.nnambu,)
@inline ndimshape(::Type{DFock}) = 1
@inline DID(i::CartesianIndex, vs::DFock) = DID(i.I...)
@inline CartesianIndex(did::DID{Int}, vs::DFock) = CartesianIndex(did.nambu)

function Base.angle(id::OID{<:Index{<:AbstractPID, DID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoord, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.iid.nambu == 1) ? phase : -phase
end
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:DID}; kwargs...) = index.pid.site
@inline script(::Val{:nambu}, index::Index{<:AbstractPID, <:DID}; kwargs...) = index.iid.nambu==2 ? "\\dagger" : ""

latexname(::Type{<:OID{<:Index{<:AbstractPID, <:DID}}}) = Symbol("OID{Index{AbstractPID, DID}}")
latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:site,)}('d'))

@testset "IID" begin
    did = DID(1)
    @test statistics(did) == statistics(typeof(did)) == :f

    did₁, did₂ = DID(1), DID(2)
    ciid = CompositeIID(did₁, did₂)
    @test length(ciid) == length(typeof(ciid)) == 2
    @test rank(ciid) == rank(typeof(ciid)) == 2
    @test iidtype(ciid, 1) == iidtype(typeof(ciid), 1) == DID{Int}
    @test iidtype(ciid, 2) == iidtype(typeof(ciid), 2) == DID{Int}
    @test ciid[1]==did₁ && ciid[2]==did₂
    @test ciid.nambus == (1, 2)
    @test string(ciid) == "DID(1) ⊗ DID(2)"
    @test did₁⊗did₂ == ciid
    @test did₁⊗ciid == CompositeIID(did₁, did₁, did₂)
    @test ciid⊗did₁ == CompositeIID(did₁, did₂, did₁)
    @test ciid⊗ciid == CompositeIID(did₁, did₂, did₁, did₂)
end

@testset "SimpleInternal" begin
    it = DFock(2)
    @test it|>eltype == DID{Int}
    @test it|>typeof|>eltype == DID{Int}
    @test it == deepcopy(it)
    @test isequal(it, deepcopy(it))
    @test it|>string == "DFock(nnambu=2)"
    @test it|>collect == [DID(1), DID(2)]
    @test statistics(it) == statistics(typeof(it)) == :f
    @test match(DID(1), it) && match(DID, DFock)
    @test filter(DID(1), it) == filter(DID, it) == it
    @test filter(DID(1), DFock) == filter(DID, DFock) == DFock
end

@testset "CompositeInternal" begin
    it₁, it₂ = DFock(2), DFock(3)
    ci = CompositeInternal(it₁, it₂)
    @test eltype(ci) == eltype(typeof(ci)) == CompositeIID{Tuple{DID{Int}, DID{Int}}}
    @test rank(ci) == rank(typeof(ci)) == 2
    @test internaltype(ci, 1) == internaltype(typeof(ci), 1) == DFock
    @test internaltype(ci, 2) == internaltype(typeof(ci), 2) == DFock
    @test ci[InternalIndex(1)]==it₁ && ci[InternalIndex(2)]==it₂
    @test string(ci) == "DFock(nnambu=2) ⊗ DFock(nnambu=3)"
    @test shape(ci) == (1:2, 1:3)
    @test ndimshape(ci) == 2
    for i = 1:2, j = 1:3
        @test CartesianIndex(CompositeIID(CartesianIndex(i, j), ci), ci) == CartesianIndex(i, j)
    end
    @test it₁⊗it₂ == ci
    @test it₁⊗ci == CompositeInternal(it₁, it₁, it₂)
    @test ci⊗it₁ == CompositeInternal(it₁, it₂, it₁)
    @test ci⊗ci == CompositeInternal(it₁, it₂, it₁, it₂)
    @test filter(DID(1), ci) == filter(DID, ci) == ci
    @test filter(DID(1), typeof(ci)) == filter(DID, typeof(ci)) == typeof(ci)
end

@testset "Hilbert" begin
    map = pid->DFock((pid.site-1)%2+1)
    hilbert = Hilbert(map, [CPID(1, 1), CPID(1, 2)])
    @test hilbert == Hilbert{DFock}(map, [CPID(1, 1), CPID(1, 2)])
    @test convert(Dict, hilbert) == Dict(CPID(1, 1)=>DFock(1), CPID(1, 2)=>DFock(2))
    reset!(hilbert, (CPID(2, 1), CPID(2, 2)))
    @test convert(Dict, hilbert) == Dict(CPID(2, 1)=>DFock(1), CPID(2, 2)=>DFock(2))

    hilbert = Hilbert(pid=>DFock(2) for pid in [PID(1), PID(2)])
    @test hilbert[PID(1)] == hilbert[PID(2)] == DFock(2)

    hilbert = Hilbert(PID(1)=>DFock(2), PID(2)=>DFock(3))
    @test hilbert[PID(1)]==DFock(2) && hilbert[PID(2)]==DFock(3)
end

@testset "Index" begin
    index = Index(CPID('S', 4), DID(1))
    @test index|>pidtype == CPID{Char}
    @test index|>typeof|>pidtype == CPID{Char}
    @test index|>iidtype == DID{Int}
    @test index|>typeof|>iidtype == DID{Int}
    @test index|>adjoint == Index(CPID('S', 4), DID(2))
    @test statistics(index) == statistics(typeof(index)) == :f
    @test isHermitian(ID(index', index)) == true
    @test isHermitian(ID(index, index)) == false
end

@testset "OID" begin
    @test contentnames(AbstractCompositeOID) == (:index,)
    @test contentnames(OID) == (:index, :rcoord, :icoord)

    oid = OID(Index(PID(1), DID(1)), [0.0, -0.0], [0.0, 0.0])
    @test indextype(oid) == indextype(typeof(oid)) == Index{PID, DID{Int}}
    @test statistics(oid) == statistics(typeof(oid)) == :f
    @test oid' == OID(Index(PID(1), DID(2)), rcoord=SVector(0.0, 0.0), icoord=SVector(0.0, 0.0))
    @test hash(oid, UInt(1)) == hash(OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 1.0)), UInt(1))
    @test propertynames(ID(oid)) == (:indexes, :rcoords, :icoords)
    @test fieldnames(OID) == (:index, :rcoord, :icoord)
    @test string(oid) == "OID(Index(PID(1), DID(1)), [0.0, 0.0], [0.0, 0.0])"
    @test ID(oid', oid)' == ID(oid', oid)
    @test isHermitian(ID(oid', oid)) == true
    @test isHermitian(ID(oid, oid)) == false
    @test oidtype(DFock, Point{2, PID, Float}, Val(:default)) == OID{Index{PID, DID{Int}}, SVector{2, Float}}
end

@testset "Operator" begin
    opt = Operator(1.0im, ID(
        OID(Index(PID(2), DID(2)), SVector(1.0, 0.0), SVector(2.0, 0.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        ))
    @test opt' == Operator(-1.0im, ID(
        OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
        OID(Index(PID(2), DID(1)), SVector(1.0, 0.0), SVector(2.0, 0.0))
        ))
    @test isHermitian(opt) == false
    @test string(opt) == "Operator(1.0im, ID(OID(Index(PID(2), DID(2)), [1.0, 0.0], [2.0, 0.0]), OID(Index(PID(1), DID(1)), [0.0, 0.0], [0.0, 0.0])))"

    opt = Operator(1.0, ID(
        OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(1.0, 1.0)),
        OID(Index(PID(1), DID(1)), SVector(0.5, 0.5), SVector(1.0, 1.0))
        ))
    @test opt' == opt
    @test isHermitian(opt) == true

    opt = Operator(1.0, ID(
        OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(1.0, 1.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.5), SVector(0.0, 1.0))
        ))
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)

    opt = Operator(1.0, ID(OID(Index(PID(1), DID(2)), SVector(0.5, 0.0), SVector(1.0, 0.0))))
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)
end

@testset "Operators" begin
    opt1 = Operator(1.0im, ID(
        OID(Index(PID(2), DID(2)), SVector(1.0, 0.0), SVector(0.0, 0.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(2.0, 0.0))
        ))
    opt2 = Operator(1.0, ID(
        OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        ))
    opts = Operators(opt1, opt2)
    @test summary(opts) == "Operators{$(opts|>valtype)}"
    @test opts' == Operators(opt1', opt2')
    @test opts'+opts == Operators(opt1, opt1', opt2*2)
    @test isHermitian(opts) == false
    @test isHermitian(opts'+opts) == true
end

@testset "Metric" begin
    valtype(Metric, AbstractOID) = AbstractOID

    index = Index(CPID('S', 4), DID(1))
    oid = OID(index, SVector(0.5, 0.0), SVector(1.0, 0.0))

    m = OIDToTuple((:scope, :site, :nambu))
    @test m == OIDToTuple(:scope, :site, :nambu)
    @test isequal(m, OIDToTuple(:scope, :site, :nambu))
    @test keys(m) == keys(typeof(m)) == (:scope, :site, :nambu)
    @test filter(≠(:nambu), m) == OIDToTuple(:scope, :site)
    @test OIDToTuple(Index{PID, DID{Int}}) == OIDToTuple(:site, :nambu)
    @test OIDToTuple(OID{Index{CPID{Int}, DID{Int}}}) == OIDToTuple(:scope, :site, :nambu)
    @test OIDToTuple(Hilbert{DFock, CPID{Int}}) == OIDToTuple(:scope, :site, :nambu)
    @test m(index) == ('S', 4, 1) == m(oid)
    @test filter(≠(:scope), m)(index) == (4, 1)
    @test filter(≠(:nambu), m)(index) == ('S', 4)
    @test filter(∉((:site, :nambu)), m)(index) == ('S',)
end

@testset "Table" begin
    @test contentnames(Table) == (:by, :contents)

    by = filter(≠(:nambu), OIDToTuple(Index{PID, DID{Int}}))

    table = Table([Index(PID(1), DID(1)), Index(PID(1), DID(2))], by)
    @test empty(table) == Table{Index{PID, DID{Int}}}(by)
    @test table[Index(PID(1), DID(1))]==1 && table[Index(PID(1), DID(2))]==1

    hilbert = Hilbert(pid=>DFock(2) for pid in [PID(1), PID(2)])
    inds1 = (Index(PID(1), iid) for iid in DFock(2))|>collect
    inds2 = (Index(PID(2), iid) for iid in DFock(2))|>collect
    @test Table(hilbert, by) == Table([inds1; inds2], by)
    @test Table(hilbert, by) == union(Table(inds1, by), Table(inds2, by))

    opt = Operator(1.0im, ID(
        OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        ))
    @test sequence(opt, table) == (1, 1)
    @test haskey(table, opt.id) == (true, true)

    table = Table(hilbert)
    @test table == Table([inds1; inds2])
    @test reset!(empty(table), [inds1; inds2]) == table
    @test reset!(empty(table), hilbert) == table
end

@testset "LaTeX" begin
    @test latexname(OID) == :OID
    @test latexformat(OID{<:Index{<:AbstractPID, <:DID}}) == LaTeX{(:nambu,), (:site,)}('d')

    latex = LaTeX{(:nambu,), (:site,)}('c', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0)))
    @test superscript(latex|>typeof) == (:nambu,)
    @test subscript(latex|>typeof) == (:site,)

    oid = OID(Index(CPID('d', 1), DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0))
    @test script(Val(:rcoord), oid) == "[0.0, 0.0]"
    @test script(Val(:icoord), oid) == "[1.0, 0.0]"
    @test script(Val(:integralicoord), oid; vectors=get(latex.options, :vectors, nothing)) == "[1, 0]"
    @test script(Val(:BD), oid, latex) == 'c'
    @test script(Val(:SP), oid, latex) == ("\\dagger",)
    @test script(Val(:SB), oid, latex) == (1,)
    @test repr(oid, latex) == "c^{\\dagger}_{1}"

    opt = Operator(1.0, ID(
        OID(Index(CPID('d', 2), DID(2)), SVector(1.0, 0.0), SVector(2.0, 0.0)),
        OID(Index(CPID('d', 2), DID(2)), SVector(1.0, 0.0), SVector(2.0, 0.0))
        ))
    @test repr(opt) == "(d^{\\dagger}_{2})^2"

    opt = Operator(1.0, ID(
        OID(Index(CPID('d', 2), DID(2)), SVector(1.0, 0.0), SVector(2.0, 0.0)),
        OID(Index(CPID('d', 1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        ))
    @test repr(opt) == "d^{\\dagger}_{2}d^{}_{1}"

    latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:site,)}('c'))
    @test repr(opt) == "c^{\\dagger}_{2}c^{}_{1}"

    latexformat(OID{<:Index{<:AbstractPID, <:DID}},  LaTeX{(:nambu,), (:site,)}('d'))
    @test repr(opt) == "d^{\\dagger}_{2}d^{}_{1}"

    latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:rcoord,)}('d'))
    @test repr(opt) == "d^{\\dagger}_{[1.0, 0.0]}d^{}_{[0.0, 0.0]}"

    latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:integralicoord,)}('d', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0))))
    @test repr(opt) == "d^{\\dagger}_{[2, 0]}d^{}_{[0, 0]}"

    opts = Operators(
            Operator(1.0-1.0im, ID(
                OID(Index(PID(2), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                )),
            Operator(-1.0, ID(
                OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                ))
            )
    latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:site,)}('c'))
    @test repr(opts) == "-c^{\\dagger}_{1}c^{}_{1}+(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}"

    latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:site,)}('d'))
    io = IOBuffer()
    show(io, MIME"text/latex"(), opt)
    @test String(take!(io)) == "\$d^{\\dagger}_{2}d^{}_{1}\$"
    show(io, MIME"text/latex"(), opts)
    @test String(take!(io)) == "\$-d^{\\dagger}_{1}d^{}_{1}+(1.0-1.0im)d^{\\dagger}_{2}d^{}_{1}\$"
end

@testset "Boundary" begin
    opt = Operator(4.5, ID(
        OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
        OID(Index(PID(2), DID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
        ))
    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    @test bound==deepcopy(bound) && isequal(bound, deepcopy(bound))
    @test angle(bound, opt) ≈ angle(opt.id, bound.vectors, bound.values) ≈ 0.6pi
    @test bound(opt) ≈ twist(opt, bound.vectors, bound.values) ≈ replace(opt, 4.5*exp(2im*pi*0.3))
    update!(bound, θ₁=0.3)
    @test bound(opt) ≈ replace(opt, 4.5*exp(2im*pi*0.5))

    @test angle(plain, opt) == 0
    @test plain(opt) == opt
    @test update!(plain) == plain
end
