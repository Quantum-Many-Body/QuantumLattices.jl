using Test
using StaticArrays: SVector
using LinearAlgebra: dot
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.Spatials: PID, Point, pidtype, rcoord, icoord
using QuantumLattices.Interfaces: decompose
using QuantumLattices.Essentials: update!, reset!
using QuantumLattices.Prerequisites.Traits: contentnames
using QuantumLattices.Mathematics.AlgebraOverFields: ID, sequence
import QuantumLattices.Essentials.DegreesOfFreedom: latexname, script

struct DID <: IID
    nambu::Int
end
@inline Base.adjoint(sl::DID) = DID(3-sl.nambu)

struct DFock <: Internal{DID}
    atom::Int
    nnambu::Int
end
@inline Base.Dims(f::DFock) = (f.nnambu,)
@inline DID(i::CartesianIndex, vs::DFock) = DID(i.I...)
@inline CartesianIndex(did::DID, vs::DFock) = CartesianIndex(did.nambu)

struct DIndex{S} <: Index{PID{S}, DID}
    scope::S
    site::Int
    nambu::Int
end
@inline Base.union(::Type{P}, ::Type{I}) where {P<:PID, I<:DID} = DIndex{fieldtype(P, :scope)}
function Base.angle(id::OID{<:DIndex}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoord, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.nambu == 1) ? phase : -phase
end
@inline script(::Val{:site}, index::DIndex; kwargs...) = index.site
@inline script(::Val{:nambu}, index::DIndex; kwargs...) = index.nambu==2 ? "\\dagger" : ""

struct DOperator{V<:Number, I<:ID{OID}} <: Operator{V, I}
    value::V
    id::I
end
latexname(::Type{<:OID{<:DIndex}}) = Symbol("OID{DIndex}")
latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:site,)}('d'))

struct IdentityMetric <: Metric end

@testset "Internal" begin
    it = DFock(1, 2)
    @test it|>eltype == DID
    @test it|>typeof|>eltype == DID
    @test it == deepcopy(it)
    @test isequal(it, deepcopy(it))
    @test it|>string == "DFock(atom=1, nnambu=2)"
    @test it|>collect == [DID(1), DID(2)]
end

@testset "Config" begin
    config = Config{DFock}(pid->DFock((pid.site-1)%2+1, 2), [PID(1, 1), PID(1, 2)])
    @test convert(Dict, config) == Dict(PID(1, 1)=>DFock(1, 2), PID(1, 2)=>DFock(2, 2))
    reset!(config, (PID(2, 1), PID(2, 2)))
    @test convert(Dict, config) == Dict(PID(2, 1)=>DFock(1, 2), PID(2, 2)=>DFock(2, 2))
end

@testset "Index" begin
    index = DIndex(PID('S', 4), DID(1))
    @test index|>pidtype == PID{Char}
    @test index|>typeof|>pidtype == PID{Char}
    @test index|>iidtype == DID
    @test index|>typeof|>iidtype == DID
    @test index|>pid == PID('S', 4)
    @test index|>iid == DID(1)
    @test index|>adjoint == DIndex('S', 4, 2)
    @test isHermitian(ID(index', index)) == true
    @test isHermitian(ID(index, index)) == false
    @test union(PID, IID) == Index{PID, IID}
    @test union(index|>pidtype, index|>iidtype) == index|>typeof
end

@testset "OID" begin
    @test contentnames(AbstractCompositeOID) == (:index,)
    @test contentnames(OID) == (:index, :rcoord, :icoord)

    oid = OID(DIndex(1, 1, 1), [0.0, -0.0], [0.0, 0.0])
    @test indextype(oid) == indextype(typeof(oid)) == DIndex{Int}
    @test oid' == OID(DIndex(1, 1, 2), rcoord=SVector(0.0, 0.0), icoord=SVector(0.0, 0.0))
    @test hash(oid, UInt(1)) == hash(OID(DIndex(1, 1, 1), SVector(0.0, 0.0), SVector(0.0, 1.0)), UInt(1))
    @test propertynames(ID(oid)) == (:indexes, :rcoords, :icoords)
    @test fieldnames(OID) == (:index, :rcoord, :icoord)
    @test string(oid) == "OID(DIndex(1, 1, 1), [0.0, 0.0], [0.0, 0.0])"
    @test ID(oid', oid)' == ID(oid', oid)
    @test isHermitian(ID(oid', oid)) == true
    @test isHermitian(ID(oid, oid)) == false
    @test oidtype(DFock, Point{2, PID{Char}}, Val(:default)) == OID{DIndex{Char}, SVector{2, Float}}
end

@testset "Operator" begin
    opt = DOperator(1.0im, ID(OID, (DIndex(1, 2, 2), DIndex(1, 1, 1)), (SVector(1.0, 0.0), SVector(0.0, 0.0)), (SVector(2.0, 0.0), SVector(0.0, 0.0))))
    @test opt' == DOperator(-1.0im, ID(OID, (DIndex(1, 1, 2), DIndex(1, 2, 1)), (SVector(0.0, 0.0), SVector(1.0, 0.0)), (SVector(0.0, 0.0), SVector(2.0, 0.0))))
    @test isHermitian(opt) == false
    @test string(opt) == "DOperator(1.0im, ID(OID(DIndex(1, 2, 2), [1.0, 0.0], [2.0, 0.0]), OID(DIndex(1, 1, 1), [0.0, 0.0], [0.0, 0.0])))"

    opt = DOperator(1.0, ID(OID, (DIndex(1, 1, 2), DIndex(1, 1, 1)), (SVector(0.5, 0.5), SVector(0.5, 0.5)), (SVector(1.0, 1.0), SVector(1.0, 1.0))))
    @test opt' == opt
    @test isHermitian(opt) == true

    opt = DOperator(1.0, ID(OID, (DIndex(1, 1, 2), DIndex(1, 1, 1)), (SVector(0.5, 0.5), SVector(0.0, 0.5)), (SVector(1.0, 1.0), SVector(0.0, 1.0))))
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)

    opt = DOperator(1.0, ID(OID(DIndex(1, 1, 2), SVector(0.5, 0.0), SVector(1.0, 0.0))))
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)
end

@testset "Operators" begin
    opt1 = DOperator(1.0im, ID(OID, (DIndex(1, 2, 2), DIndex(1, 1, 1)), (SVector(1.0, 0.0), SVector(0.0, 0.0)), (SVector(2.0, 0.0), SVector(0.0, 0.0))))
    opt2 = DOperator(1.0, ID(OID, (DIndex(1, 1, 2), DIndex(1, 1, 1)), (SVector(0.0, 0.0), SVector(0.0, 0.0)), (SVector(0.0, 0.0), SVector(0.0, 0.0))))
    opts = Operators(opt1, opt2)
    @test summary(opts) == "Operators{$(opts|>valtype)}"
    @test opts' == Operators(opt1', opt2')
    @test opts'+opts == Operators(opt1, opt1', opt2*2)
    @test isHermitian(opts) == false
    @test isHermitian(opts'+opts) == true
end

@testset "Metric" begin
    valtype(Metric, AbstractOID) = AbstractOID

    m = OIDToTuple((:scope, :site, :nambu))
    @test m == OIDToTuple(:scope, :site, :nambu)
    @test isequal(m, OIDToTuple(:scope, :site, :nambu))
    @test keys(m) == keys(typeof(m)) == (:scope, :site, :nambu)
    @test filter(≠(:nambu), m) == OIDToTuple(:scope, :site)

    @test OIDToTuple(DIndex) == OIDToTuple(:scope, :site, :nambu)
    @test OIDToTuple(OID{DIndex}) == OIDToTuple(:scope, :site, :nambu)
    @test OIDToTuple(Config{DFock, PID{Int}}) == OIDToTuple(:scope, :site, :nambu)

    n = IdentityMetric()
    index = DIndex(PID('S', 4), DID(1))
    oid = OID(index, SVector(0.5, 0.0), SVector(1.0, 0.0))

    @test n(index) == index && n(oid) == oid
    @test m(index) == ('S', 4, 1) == m(oid)

    @test filter(≠(:scope), m)(index) == (4, 1)
    @test filter(≠(:nambu), m)(index) == ('S', 4)
    @test filter(∉((:site, :nambu)), m)(index) == ('S',)
end

@testset "Table" begin
    @test contentnames(Table) == (:by, :contents)

    by = filter(≠(:nambu), OIDToTuple(DIndex))

    table = Table([DIndex(1, 1, 1), DIndex(1, 1, 2)], by)
    @test empty(table) == Table{DIndex}(by)
    @test table[DIndex(1, 1, 1)]==1 && table[DIndex(1, 1, 2)]==1

    config = Config{DFock}(pid->DFock((pid.site-1)%2+1, 2), [PID(1, 1), PID(1, 2)])
    inds1 = (DIndex(PID(1, 1), iid) for iid in DFock(1, 2))|>collect
    inds2 = (DIndex(PID(1, 2), iid) for iid in DFock(2, 2))|>collect
    @test Table(config, by) == Table([inds1; inds2], by)
    @test Table(config, by) == union(Table(inds1, by), Table(inds2, by))

    opt = DOperator(1.0im, ID(OID, (DIndex(1, 1, 2), DIndex(1, 1, 1)), (SVector(0.0, 0.0), SVector(0.0, 0.0)), (SVector(1.0, 0.0), SVector(0.0, 0.0))))
    @test sequence(opt, table) == (1, 1)
    @test haskey(table, opt.id) == (true, true)

    table = Table(config)
    @test table == Table([inds1; inds2])
    @test reset!(empty(table), [inds1; inds2]) == table
    @test reset!(empty(table), config) == table
end

@testset "LaTeX" begin
    @test latexname(OID) == :OID
    @test latexformat(OID{<:DIndex}) == LaTeX{(:nambu,), (:site,)}('d')

    latex = LaTeX{(:nambu,), (:site,)}('c', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0)))
    @test superscript(latex|>typeof) == (:nambu,)
    @test subscript(latex|>typeof) == (:site,)

    oid = OID(DIndex('d', 1, 2), SVector(0.0, 0.0), SVector(1.0, 0.0))
    @test script(Val(:rcoord), oid) == "[0.0, 0.0]"
    @test script(Val(:icoord), oid) == "[1.0, 0.0]"
    @test script(Val(:integeralicoord), oid; vectors=get(latex.options, :vectors, nothing)) == "[1, 0]"
    @test script(Val(:BD), oid, latex) == 'c'
    @test script(Val(:SP), oid, latex) == ("\\dagger",)
    @test script(Val(:SB), oid, latex) == (1,)
    @test repr(oid, latex) == "c^{\\dagger}_{1}"

    opt = DOperator(1.0, ID(OID, (DIndex('d', 2, 2), DIndex('d', 1, 1)), (SVector(1.0, 0.0), SVector(0.0, 0.0)), (SVector(2.0, 0.0), SVector(0.0, 0.0))))
    @test repr(opt) == "d^{\\dagger}_{2}d^{}_{1}"

    latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:site,)}('c'))
    @test repr(opt) == "c^{\\dagger}_{2}c^{}_{1}"

    latexformat(OID{<:DIndex},  LaTeX{(:nambu,), (:site,)}('d'))
    @test repr(opt) == "d^{\\dagger}_{2}d^{}_{1}"

    latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:rcoord,)}('d'))
    @test repr(opt) == "d^{\\dagger}_{[1.0, 0.0]}d^{}_{[0.0, 0.0]}"

    latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:integeralicoord,)}('d', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0))))
    @test repr(opt) == "d^{\\dagger}_{[2, 0]}d^{}_{[0, 0]}"

    opts = Operators(
            DOperator(1.0-1.0im, ID(OID, (DIndex('d', 2, 2), DIndex('d', 1, 1)), (SVector(0.0, 0.0), SVector(0.0, 0.0)), (SVector(0.0, 0.0), SVector(0.0, 0.0)))),
            DOperator(-1.0, ID(OID, (DIndex('d', 1, 2), DIndex('d', 1, 1)), (SVector(0.0, 0.0), SVector(0.0, 0.0)), (SVector(0.0, 0.0), SVector(0.0, 0.0))))
            )
    latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:site,)}('c'))
    @test repr(opts) == "-c^{\\dagger}_{1}c^{}_{1}+(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}"

    latexformat(OID{<:DIndex}, LaTeX{(:nambu,), (:site,)}('d'))
    io = IOBuffer()
    show(io, MIME"text/latex"(), opt)
    @test String(take!(io)) == "\$d^{\\dagger}_{2}d^{}_{1}\$"
    show(io, MIME"text/latex"(), opts)
    @test String(take!(io)) == "\$-d^{\\dagger}_{1}d^{}_{1}+(1.0-1.0im)d^{\\dagger}_{2}d^{}_{1}\$"
end

@testset "Boundary" begin
    opt = DOperator(4.5, ID(OID, (DIndex('a', 1, 2), DIndex('b', 2, 1)), (SVector(0.5, 0.5), SVector(1.5, 1.5)), (SVector(0.0, 0.0), SVector(1.0, 1.0))))
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
