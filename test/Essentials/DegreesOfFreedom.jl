using Test
using Printf: @printf, @sprintf
using StaticArrays: SVector
using LinearAlgebra: dot
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.QuantumOperators: ID, Operator, Operators, sequence, id, idtype, LaTeX, latexformat
using QuantumLattices.Essentials.Spatials: AbstractPID, PID, CPID, Point, Bond, Bonds, Lattice, pidtype, rcoord, icoord
using QuantumLattices.Essentials: kind, update!, reset!
using QuantumLattices.Interfaces: decompose, rank, expand, ⊗
using QuantumLattices.Prerequisites: Float, decimaltostr
using QuantumLattices.Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
using QuantumLattices.Prerequisites.CompositeStructures: NamedContainer

import QuantumLattices.Essentials.QuantumOperators: latexname, script
import QuantumLattices.Essentials.DegreesOfFreedom: ishermitian, statistics, couplingcenters, abbr
import QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

struct DID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.show(io::IO, did::DID) = @printf io "DID(%s)" did.nambu
@inline Base.adjoint(sl::DID{Int}) = DID(3-sl.nambu)
@inline statistics(::Type{<:DID}) = :f

struct DFock <: SimpleInternal{DID{Int}}
    nnambu::Int
end
@inline shape(f::DFock) = (1:f.nnambu,)
@inline ndimshape(::Type{DFock}) = 1
@inline DID(i::CartesianIndex, vs::DFock) = DID(i.I...)
@inline CartesianIndex(did::DID{Int}, vs::DFock) = CartesianIndex(did.nambu)
@inline shape(iidspace::IIDSpace{DID{Symbol}, DFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{DID{Int}, DFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)
@inline shape(iidspace::IIDSpace{DID{Symbol}, DFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{DID{Int}, DFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

function Base.angle(id::OID{<:Index{<:AbstractPID, DID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoord, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.iid.nambu == 1) ? phase : -phase
end
@inline script(::Val{:site}, index::Index{<:AbstractPID, <:DID}; kwargs...) = index.pid.site
@inline script(::Val{:nambu}, index::Index{<:AbstractPID, <:DID}; kwargs...) = index.iid.nambu==2 ? "\\dagger" : ""

@inline latexname(::Type{<:OID{<:Index{<:AbstractPID, <:DID}}}) = Symbol("OID{Index{AbstractPID, DID}}")
latexformat(OID{<:Index{<:AbstractPID, <:DID}}, LaTeX{(:nambu,), (:site,)}('d'))

const DCoupling{V, I<:ID{DID}, C<:Subscripts, CI<:SubscriptsID} = Coupling{V, I, C, CI}
@inline Base.repr(tc::(Coupling{V, <:ID{DID}} where V)) = @sprintf "%s ph(%s)" decimaltostr(tc.value) join(tc.cid.nambus, ", ")
@inline couplingcenters(::(Coupling{V, <:ID{DID}} where V), ::Bond, ::Val) = (1, 2)
@inline DCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(DID, nambus), Subscripts((nambu=Subscript(nambus),)))
@inline DCoupling(value, nambus::Subscript) = Coupling(value, ID(DID, convert(Tuple, nambus)), Subscripts((nambu=nambus,)))

@inline abbr(::Type{<:Term{:Mu}}) = :mu
@inline abbr(::Type{<:Term{:Hp}}) = :hp
@inline ishermitian(::Type{<:Term{:Mu}}) = true
@inline ishermitian(::Type{<:Term{:Hp}}) = false

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

@testset "Index" begin
    @test parameternames(Index) == (:pid, :iid)
    @test isparameterbound(Index, Val(:pid), PID) == false
    @test isparameterbound(Index, Val(:iid), DID) == true

    index = Index(CPID('S', 4), DID(1))
    @test index|>pidtype == CPID{Char}
    @test index|>typeof|>pidtype == CPID{Char}
    @test index|>iidtype == DID{Int}
    @test index|>typeof|>iidtype == DID{Int}
    @test index|>adjoint == Index(CPID('S', 4), DID(2))
    @test statistics(index) == statistics(typeof(index)) == :f
    @test ishermitian(ID(index', index)) == true
    @test ishermitian(ID(index, index)) == false
end

@testset "OID" begin
    @test contentnames(CompositeOID) == (:index,)
    @test parameternames(CompositeOID) == (:index,)
    @test isparameterbound(CompositeOID, Val(:index), Index) == true
    @test isparameterbound(CompositeOID, Val(:index), Index{PID, DID{Int}}) == false

    @test contentnames(OID) == (:index, :rcoord, :icoord)
    @test parameternames(OID) == (:index, :coord)
    @test isparameterbound(OID, Val(:index), Index) == true
    @test isparameterbound(OID, Val(:index), Index{PID, DID{Int}}) == false
    @test isparameterbound(OID, Val(:coord), SVector) == true
    @test isparameterbound(OID, Val(:coord), SVector{2, Float}) == false

    oid = OID(Index(PID(1), DID(1)), [0.0, -0.0], [0.0, 0.0])
    @test indextype(oid) == indextype(typeof(oid)) == Index{PID, DID{Int}}
    @test statistics(oid) == statistics(typeof(oid)) == :f
    @test oid' == OID(Index(PID(1), DID(2)), rcoord=SVector(0.0, 0.0), icoord=SVector(0.0, 0.0))
    @test hash(oid, UInt(1)) == hash(OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 1.0)), UInt(1))
    @test propertynames(ID(oid)) == (:indexes, :rcoords, :icoords)
    @test fieldnames(OID) == (:index, :rcoord, :icoord)
    @test string(oid) == "OID(Index(PID(1), DID(1)), [0.0, 0.0], [0.0, 0.0])"
    @test ID(oid', oid)' == ID(oid', oid)
    @test ishermitian(ID(oid', oid)) == true
    @test ishermitian(ID(oid, oid)) == false
    @test oidtype(DFock, Point{2, PID, Float}, Val(:default)) == OID{Index{PID, DID{Int}}, SVector{2, Float}}
end

@testset "Operator" begin
    opt = Operator(1.0,
        OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(1.0, 1.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.5), SVector(0.0, 1.0))
        )
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)

    opt = Operator(1.0, OID(Index(PID(1), DID(2)), SVector(0.5, 0.0), SVector(1.0, 0.0)))
    @test rcoord(opt) == SVector(0.5, 0.0)
    @test icoord(opt) == SVector(1.0, 0.0)
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

@testset "IIDSpace" begin
    DID₁, DID₂, it = DID(2), DID(:σ), DFock(2)
    iidspace = IIDSpace(DID₁⊗DID₂, it⊗it)
    @test eltype(iidspace) == eltype(typeof(iidspace)) == CompositeIID{Tuple{DID{Int}, DID{Int}}}
    @test kind(iidspace) == kind(typeof(iidspace)) == :info
    @test shape(iidspace) == (2:2, 1:2)
    @test ndimshape(iidspace) == ndimshape(typeof(iidspace)) == 2
    for i = 1:length(iidspace)
        iid = iidspace[i]
        @test iidspace[CartesianIndex(iid, iidspace)] == iid
    end
    @test expand((DID₁, DID₂), (it, it)) == iidspace
    @test collect(iidspace) == [DID(2)⊗DID(1), DID(2)⊗DID(2)]
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

    opt = Operator(1.0im,
        OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0)),
        OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        )
    @test sequence(opt, table) == (1, 1)
    @test haskey(table, opt.id) == (true, true)

    table = Table(hilbert)
    @test table == Table([inds1; inds2])
    @test reset!(empty(table), [inds1; inds2]) == table
    @test reset!(empty(table), hilbert) == table
end

@testset "Subscript" begin
    sub = Subscript(4)
    @test contentnames(typeof(sub)) == (:contents, :rep, :constraint)
    @test getcontent(sub, :contents) == sub.pattern
    @test sub==deepcopy(sub) && isequal(sub, deepcopy(sub))
    @test string(sub) == "[* * * *]"
    @test rank(sub) == rank(typeof(sub)) == 4
    @test match(sub, (2, 2, 2, 2))

    sub = Subscript((1, 2, 3, 4))
    @test string(sub) == "[1 2 3 4]"
    @test sub == subscript"[1 2 3 4]" == subscript"[1 2 3 4]; false"
    @test rank(sub) == 4
    @test match(sub, (1, 3, 3, 4))

    sub = Subscript((1, 2, 2, 1), true)
    @test string(sub) == "[1 2 2 1]"
    @test sub == subscript"[1 2 2 1]; true"
    @test rank(sub) == 4
    @test match(sub, (1, 2, 2, 1)) && !match(sub, (1, 1, 2, 1))

    sub = subscript"[x₁ x₂ x₁ x₂]"
    @test string(sub) == "[x₁ x₂ x₁ x₂]"
    @test rank(sub) == 4
    @test match(sub, (1, 2, 1, 2)) && !match(sub, (1, 2, 2, 1))

    sub = subscript"[x₁ 4 4 x₂](x₁ < x₂)"
    @test string(sub) == "[x₁ 4 4 x₂](x₁ < x₂)"
    @test rank(sub) == rank(typeof(sub)) == 4
    @test match(sub, (1, 4, 4, 2)) && !match(sub, (2, 4, 4, 1))
end

@testset "Subscripts" begin
    subscripts = Subscripts((nambu=subscript"[x y](x > y)",), (nambu=subscript"[x x y y](x < y)",))
    @test string(subscripts) == "nambu[x y](x > y) × nambu[x x y y](x < y)"
    @test repr(subscripts, 1:2, :nambu) == "[x y](x > y)×[x x y y](x < y)"
    @test repr(subscripts, 1, :nambu) == "[x y](x > y)"
    @test repr(subscripts, 2, :nambu) == "[x x y y](x < y)"

    @test rank(subscripts) == rank(typeof(subscripts)) == 6
    @test rank(subscripts, 1) == rank(typeof(subscripts), 1) == 2
    @test rank(subscripts, 2) == rank(typeof(subscripts), 2) == 4
    @test match(subscripts, (DID(2), DID(1), DID(2), DID(2), DID(3), DID(3)))
    @test match(subscripts, DID(2)⊗DID(1)⊗DID(2)⊗DID(2)⊗DID(3)⊗DID(3))
    @test !match(subscripts, DID(1)⊗DID(2)⊗DID(2)⊗DID(2)⊗DID(3)⊗DID(3))
    @test !match(subscripts, DID(2)⊗DID(1)⊗DID(3)⊗DID(3)⊗DID(2)⊗DID(2))
    @test subscripts == Subscripts((nambu=subscript"[x y](x > y)",))*Subscripts((nambu=subscript"[x x y y](x < y)",))

    subscriptsid = SubscriptsID(subscripts)
    @test subscriptsid == SubscriptsID((1:2=>("[x y](x > y)",), 3:6=>("[x x y y](x < y)",)))
    @test idtype(subscripts) == idtype(typeof(subscripts)) == SubscriptsID{NTuple{2, Pair{UnitRange{Int}, Tuple{String}}}}
end

@testset "couplings" begin
    @test parameternames(Coupling) == (:value, :cid, :subscripts, :subscriptsid)
    @test contentnames(Coupling) == (:value, :id, :subscripts)
    @test isparameterbound(Coupling, :cid, Tuple{DID{Int}}) == false
    @test isparameterbound(Coupling, :cid, ID{DID{Int}}) == true
    @test isparameterbound(Coupling, :subscripts, Subscripts{Tuple{NamedTuple{(:nambu,), Tuple{Subscript{Tuple{Int}, typeof(noconstrain)}}}}}) == false
    @test isparameterbound(Coupling, :subscripts, Subscripts) == true
    @test isparameterbound(Coupling, :subscriptsid, SubscriptsID{Tuple{}}) == false
    @test isparameterbound(Coupling, :subscriptsid, SubscriptsID) == true

    tc = DCoupling(2.0, (2,))
    @test id(tc) == ID(CompositeIID(tc.cid), SubscriptsID(tc.subscripts))
    @test rank(tc) == rank(typeof(tc)) == 1
    @test tc == Coupling(2.0, id(tc), tc.subscripts)
    @test ID{SimpleIID}(tc) == tc.cid
    @test Subscripts(tc) == tc.subscripts

    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>DFock(2))
    @test couplingcenters(tc, point, Val(:Mu)) == (1,)
    @test couplingpoints(tc, point, Val(:Mu)) == (point,)
    @test couplinginternals(tc, point, hilbert, Val(:Mu)) == (DFock(2),)

    tc₁ = DCoupling(1.5, (1, 2))
    tc₂ = DCoupling(2.0, subscript"[a b](a < b)")
    @test tc₁*tc₂ == Coupling(3.0, ID(DID(1), DID(2), DID(:a), DID(:b)), Subscripts((nambu=Subscript((1, 2)),), (nambu=subscript"[a b](a < b)",)))

    ex = expand(tc₁, point, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Tuple{Float64, NTuple{2, OID{Index{CPID{Int}, DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        (1.5, ID(OID(Index(CPID(1, 1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                 OID(Index(CPID(1, 1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                 ))
        ]

    ex = expand(tc₂, point, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Tuple{Float64, NTuple{2, OID{Index{CPID{Int}, DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        (2.0, ID(OID(Index(CPID(1, 1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                    OID(Index(CPID(1, 1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                    ))
        ]

    @test @couplings(DCoupling(1.0, (1, 1))) == Couplings(DCoupling(1.0, (1, 1)))
    @test @couplings(DCoupling(1.0, (1,))+DCoupling(1.0, (2,))) == DCoupling(1.0, (1,))+DCoupling(1.0, (2,))
end

@testset "TermFunction" begin
    ta = TermAmplitude()
    @test ta() == 1

    ta = TermAmplitude(x->x+3.0)
    @test ta(1.0) == 4.0

    tcs = DCoupling(1.0, (1, 1)) + DCoupling(2.0, (2, 2))
    termcouplings = TermCouplings(tcs)
    @test termcouplings == deepcopy(TermCouplings(tcs))
    @test isequal(termcouplings, deepcopy(TermCouplings(tcs)))
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typeof(tcs)
    @test termcouplings() == tcs

    fx = i -> i%2==1 ? (DCoupling(1.0, (1, 1))+DCoupling(1.0, (2, 2))) : (DCoupling(1.0, (2, 1))+DCoupling(1.0, (1, 2)))
    termcouplings = TermCouplings(fx)
    @test termcouplings == TermCouplings{typejoin(typeof(fx(1)), typeof(fx(2)))}(fx)
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typejoin(typeof(fx(1)), typeof(fx(2)))
    @test termcouplings(1) == fx(1)
    @test termcouplings(2) == fx(2)

    @test ismodulatable(TermModulate{Val{true}}) == ismodulatable(TermModulate{<:Function}) == true
    @test ismodulatable(TermModulate{Val{false}}) == false
    termmodulate = TermModulate(:t)
    @test termmodulate(t=1) == 1
    @test isnothing(termmodulate(mu=1))
    @test isnothing(termmodulate())
    @test ismodulatable(termmodulate) == true

    termmodulate = TermModulate(:t, t->t*2.0)
    @test ismodulatable(termmodulate) == true
    @test termmodulate(2) == 4
end

@testset "Term" begin
    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>DFock(2))
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=(bond->3), modulate=false)
    @test term|>kind == term|>typeof|>kind == :Mu
    @test term|>id == term|>typeof|>id == :mu
    @test term|>valtype == term|>typeof|>valtype == Float
    @test term|>rank == term|>typeof|>rank == 2
    @test term|>abbr == term|>typeof|>abbr == :mu
    @test term|>ismodulatable == term|>typeof|>ismodulatable == false
    @test term|>ishermitian == term|>typeof|>ishermitian == true
    @test term == deepcopy(term)
    @test isequal(term, deepcopy(term))
    @test string(term) == "Mu{2}(id=mu, value=1.5, bondkind=0, factor=1.0)"
    @test repr(term, point, hilbert) == "mu: 4.5 ph(2, 1)"
    @test +term == term
    @test -term == term*(-1) == replace(term, factor=-term.factor)
    @test 2*term == term*2 == replace(term, factor=term.factor*2)
    @test term/2 == term*(1/2) == replace(term, factor=term.factor/2)

    p1 = Point(PID(1), (0.0, 0.0), (0.0, 0.0))
    p2 = Point(PID(2), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(pid=>DFock(2) for pid in [PID(1), PID(2)])
    term = Term{:Mu}(:mu, 1.5, 0,
        couplings=bond->bond.pid.site%2==0 ? Couplings(DCoupling(1.0, (2, 2))) : Couplings(DCoupling(1.0, (1, 1))),
        amplitude=bond->3,
        modulate=true
    )
    @test term|>ismodulatable == term|>typeof|>ismodulatable == true
    @test repr(term, p1, hilbert) == "mu: 4.5 ph(1, 1)"
    @test repr(term, p2, hilbert) == "mu: 4.5 ph(2, 2)"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test term.modulate(mu=4.0) == 4.0
    @test isnothing(term.modulate(t=1.0))
    @test update!(term, mu=4.25) == replace(term, value=4.25)
    @test term.value == 4.25
end

@testset "expand" begin
    point = Point(PID(1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>DFock(2))
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(+2.25,
            OID(Index(PID(1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            OID(Index(PID(1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    @test expand(term, point, hilbert, half=true) == operators
    @test expand(term, point, hilbert, half=false) == operators*2

    bond = Bond(1, Point(PID(2), (1.5, 1.5), (1.0, 1.0)), Point(PID(1), (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(pid=>DFock(2) for pid in [PID(1), PID(2)])
    term = Term{:Hp}(:t, 1.5, 1, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(4.5,
            OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
            OID(Index(PID(2), DID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
            )
        )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    lattice = Lattice(:Tuanzi, [Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))], vectors=[[1.0, 0.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert(pid=>DFock(2) for pid in lattice.pids)
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(+2.25,
            OID(Index(CPID(1, 1), DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            OID(Index(CPID(1, 1), DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    @test expand(term, bonds, hilbert, half=true) == operators
    @test expand(term, bonds, hilbert, half=false) == operators*2
end

@testset "script" begin
    oid = OID(Index(CPID('d', 1), DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0))
    latex = LaTeX{(:nambu,), (:site,)}('c', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0)))
    @test script(Val(:rcoord), oid) == "[0.0, 0.0]"
    @test script(Val(:icoord), oid) == "[1.0, 0.0]"
    @test script(Val(:integralicoord), oid; vectors=get(latex.options, :vectors, nothing)) == "[1, 0]"
    @test script(Val(:site), oid) == 1
    @test script(Val(:nambu), oid) == "\\dagger"
end

@testset "Boundary" begin
    opt = Operator(4.5,
        OID(Index(PID(1), DID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
        OID(Index(PID(2), DID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
        )
    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    @test keys(bound) == keys(typeof(bound)) == (:θ₁, :θ₂)
    @test bound==deepcopy(bound) && isequal(bound, deepcopy(bound))
    @test angle(bound, opt) ≈ angle(opt.id, bound.vectors, bound.values) ≈ 0.6pi
    @test bound(opt) ≈ twist(opt, bound.vectors, bound.values) ≈ replace(opt, 4.5*exp(2im*pi*0.3))
    update!(bound, θ₁=0.3)
    @test bound(opt) ≈ replace(opt, 4.5*exp(2im*pi*0.5))

    @test angle(plain, opt) == 0
    @test plain(opt) == opt
    @test update!(plain) == plain
end
