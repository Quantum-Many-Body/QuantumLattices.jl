using LinearAlgebra: dot
using Printf: @printf, @sprintf
using QuantumLattices.Essentials: kind, reset!, update!
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.QuantumOperators: ID, LaTeX, Operator, Operators, id, latexformat, sequence
using QuantumLattices.Essentials.Spatials: Bond, Lattice, Point, bonds, decompose, icoordinate, rcoordinate
using QuantumLattices.Interfaces: ⊕, ⊗, expand, rank
using QuantumLattices.Prerequisites: Float, decimaltostr
using QuantumLattices.Prerequisites.Traits: contentnames, getcontent, isparameterbound, parameternames, reparameter
using StaticArrays: SVector

import QuantumLattices.Essentials.DegreesOfFreedom: statistics
import QuantumLattices.Essentials.QuantumOperators: latexname, script
import QuantumLattices.Prerequisites.VectorSpaces: shape

struct DID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.show(io::IO, did::DID) = @printf io "DID(%s)" did.nambu
@inline Base.adjoint(sl::DID{Int}) = DID(3-sl.nambu)
@inline statistics(::Type{<:DID}) = :f
@inline DID(did::DID, ::CompositeInternal) = did
function Base.angle(id::CompositeIndex{<:Index{DID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoordinate, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.iid.nambu == 1) ? phase : -phase
end

struct DFock <: SimpleInternal{DID{Int}}
    nnambu::Int
end
@inline shape(f::DFock) = (1:f.nnambu,)
@inline DID(i::CartesianIndex, vs::DFock) = DID(i.I...)
@inline CartesianIndex(did::DID{Int}, vs::DFock) = CartesianIndex(did.nambu)
@inline shape(iidspace::IIDSpace{DID{Symbol}, DFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{DID{Int}, DFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)
@inline shape(iidspace::IIDSpace{DID{Symbol}, DFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{DID{Int}, DFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

@inline script(::Val{:site}, index::Index{<:DID}; kwargs...) = index.site
@inline script(::Val{:nambu}, index::Index{<:DID}; kwargs...) = index.iid.nambu==2 ? "\\dagger" : ""

@inline latexname(::Type{<:CompositeIndex{<:Index{<:DID}}}) = Symbol("CompositeIndex{Index{DID}}")
latexformat(CompositeIndex{<:Index{<:DID}}, LaTeX{(:nambu,), (:site,)}('d'))

const DCoupling{V, I<:ID{DID}, C<:Subscripts} = Coupling{V, I, C}
@inline Base.repr(tc::(Coupling{V, <:ID{DID}} where V)) = @sprintf "%s*ph[%s]" decimaltostr(tc.value) join(tc.iids.nambus, " ")
@inline DCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(DID, nambus), Subscripts((nambu=Subscript(nambus),)))
@inline DCoupling(value, nambus::Subscript) = Coupling(value, ID(DID, Tuple(nambus)), Subscripts((nambu=nambus,)))

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

    ci = CompositeInternal{:⊕}(it₁, it₂)
    @test eltype(ci) == eltype(typeof(ci)) == DID{Int}
    @test rank(ci) == rank(typeof(ci)) == 2
    @test string(ci) == "DFock(nnambu=2) ⊕ DFock(nnambu=3)"
    for i = 1:2
        @test ci[i] == it₁[i]
    end
    for i = 1:3
        @test ci[2+i] == it₂[i]
    end
    @test it₁⊕it₂ == ci
    @test it₁⊕ci == CompositeInternal{:⊕}(it₁, it₁, it₂)
    @test ci⊕it₁ == CompositeInternal{:⊕}(it₁, it₂, it₁)
    @test ci⊕ci == CompositeInternal{:⊕}(it₁, it₂, it₁, it₂)
    @test filter(DID(1), ci) == filter(DID, ci) == ci
    @test filter(DID(1), typeof(ci)) == filter(DID, typeof(ci)) == typeof(ci)

    ci = CompositeInternal{:⊗}(it₁, it₂)
    @test eltype(ci) == eltype(typeof(ci)) == CompositeIID{Tuple{DID{Int}, DID{Int}}}
    @test rank(ci) == rank(typeof(ci)) == 2
    @test string(ci) == "DFock(nnambu=2) ⊗ DFock(nnambu=3)"
    count = 1
    for i = 1:3
        for j = 1:2
            @test ci[count] == CompositeIID(it₁[j], it₂[i])
            count += 1
        end
    end
    @test it₁⊗it₂ == ci
    @test it₁⊗ci == CompositeInternal{:⊗}(it₁, it₁, it₂)
    @test ci⊗it₁ == CompositeInternal{:⊗}(it₁, it₂, it₁)
    @test ci⊗ci == CompositeInternal{:⊗}(it₁, it₂, it₁, it₂)
    @test filter(DID(1), ci) == filter(DID, ci) == ci
    @test filter(DID(1), typeof(ci)) == filter(DID, typeof(ci)) == typeof(ci)
end

@testset "Index" begin
    @test parameternames(Index) == (:iid,)
    @test isparameterbound(Index, Val(:iid), DID) == true

    index = Index(4, DID(1))
    @test index|>iidtype == DID{Int}
    @test index|>typeof|>iidtype == DID{Int}
    @test index|>adjoint == Index(4, DID(2))
    @test statistics(index) == statistics(typeof(index)) == :f
    @test ishermitian(ID(index', index)) == true
    @test ishermitian(ID(index, index)) == false
end

@testset "CompositeIndex" begin
    @test contentnames(AbstractCompositeIndex) == (:index,)
    @test parameternames(AbstractCompositeIndex) == (:index,)
    @test isparameterbound(AbstractCompositeIndex, Val(:index), Index) == true
    @test isparameterbound(AbstractCompositeIndex, Val(:index), Index{DID{Int}}) == false

    @test contentnames(CompositeIndex) == (:index, :rcoordinate, :icoordinate)
    @test parameternames(CompositeIndex) == (:index, :coordination)
    @test isparameterbound(CompositeIndex, Val(:index), Index) == true
    @test isparameterbound(CompositeIndex, Val(:index), Index{DID{Int}}) == false
    @test isparameterbound(CompositeIndex, Val(:coordination), SVector) == true
    @test isparameterbound(CompositeIndex, Val(:coordination), SVector{2, Float}) == false

    index = CompositeIndex(Index(1, DID(1)), [0.0, -0.0], [0.0, 0.0])
    @test indextype(index) == indextype(typeof(index)) == Index{DID{Int}}
    @test statistics(index) == statistics(typeof(index)) == :f
    @test index' == CompositeIndex(Index(1, DID(2)), rcoordinate=SVector(0.0, 0.0), icoordinate=SVector(0.0, 0.0))
    @test hash(index, UInt(1)) == hash(CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 1.0)), UInt(1))
    @test propertynames(ID(index)) == (:indexes, :rcoordinates, :icoordinates)
    @test fieldnames(CompositeIndex) == (:index, :rcoordinate, :icoordinate)
    @test string(index) == "CompositeIndex(Index(1, DID(1)), [0.0, 0.0], [0.0, 0.0])"
    @test ID(index', index)' == ID(index', index)
    @test ishermitian(ID(index', index)) == true
    @test ishermitian(ID(index, index)) == false
    @test indextype(DFock, Point{2, Float}, Val(:default)) == CompositeIndex{Index{DID{Int}}, SVector{2, Float}}
end

@testset "Operator" begin
    opt = Operator(1.0,
        CompositeIndex(Index(1, DID(2)), SVector(0.5, 0.5), SVector(1.0, 1.0)),
        CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.5), SVector(0.0, 1.0))
        )
    @test rcoordinate(opt) == SVector(0.5, 0.0)
    @test icoordinate(opt) == SVector(1.0, 0.0)

    opt = Operator(1.0, CompositeIndex(Index(1, DID(2)), SVector(0.5, 0.0), SVector(1.0, 0.0)))
    @test rcoordinate(opt) == SVector(0.5, 0.0)
    @test icoordinate(opt) == SVector(1.0, 0.0)
end

@testset "Hilbert" begin
    map = site->DFock((site-1)%2+1)
    hilbert = Hilbert(map, [1, 2])
    @test hilbert == Hilbert{DFock}(map, [1, 2])
    reset!(hilbert, [3, 4])
    @test hilbert == Hilbert{DFock}(map, [3, 4])

    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    @test hilbert[1] == hilbert[2] == DFock(2)

    hilbert = Hilbert(1=>DFock(2), 2=>DFock(3))
    @test hilbert[1]==DFock(2) && hilbert[2]==DFock(3)
end

@testset "IIDSpace" begin
    DID₁, DID₂, it = DID(2), DID(:σ), DFock(2)
    iidspace = IIDSpace(DID₁⊗DID₂, it⊗it)
    @test eltype(iidspace) == eltype(typeof(iidspace)) == CompositeIID{Tuple{DID{Int}, DID{Int}}}
    @test kind(iidspace) == kind(typeof(iidspace)) == :info
    @test length(iidspace) == 2
    @test iidspace[1]==DID(2)⊗DID(1) && iidspace[2]==DID(2)⊗DID(2)
    @test expand((DID₁, DID₂), (it, it)) == iidspace
    @test collect(iidspace) == [DID(2)⊗DID(1), DID(2)⊗DID(2)]
end

@testset "Metric" begin
    valtype(Metric, AbstractOID) = AbstractOID

    index = CompositeIndex(Index(4, DID(1)), SVector(0.5, 0.0), SVector(1.0, 0.0))

    m = OperatorUnitToTuple((:site, :nambu))
    @test m == OperatorUnitToTuple(:site, :nambu)
    @test isequal(m, OperatorUnitToTuple(:site, :nambu))
    @test keys(m) == keys(typeof(m)) == (:site, :nambu)
    @test OperatorUnitToTuple(Index{DID{Int}}) == OperatorUnitToTuple(:site, :nambu)
    @test OperatorUnitToTuple(Hilbert{DFock}) == OperatorUnitToTuple(:site, :nambu)
    @test m(index.index) == (4, 1) == m(index)
end

@testset "Table" begin
    @test contentnames(Table) == (:by, :contents)

    by = OperatorUnitToTuple(:site)

    table = Table([Index(1, DID(1)), Index(1, DID(2))], by)
    @test empty(table) == Table{Index{DID{Int}}}(by)
    @test table[Index(1, DID(1))]==1 && table[Index(1, DID(2))]==1

    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    inds₁ = (Index(1, iid) for iid in DFock(2))|>collect
    inds₂ = (Index(2, iid) for iid in DFock(2))|>collect
    @test Table(hilbert, by) == Table([inds₁; inds₂], by)
    @test Table(hilbert, by) == union(Table(inds₁, by), Table(inds₂, by))

    opt = Operator(1.0im,
        CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0)),
        CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
        )
    @test sequence(opt, table) == (1, 1)
    @test haskey(table, opt.id) == (true, true)

    table = Table(hilbert)
    @test table == Table([inds₁; inds₂])
    @test reset!(empty(table), [inds₁; inds₂]) == table
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
    @test subscripts.rep == (1:2=>("[x y](x > y)",), 3:6=>("[x x y y](x < y)",))
    @test hash(subscripts) == hash(subscripts.rep)
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
end

@testset "couplings" begin
    @test parameternames(Coupling) == (:value, :iids, :subscripts)
    @test isparameterbound(Coupling, :iids, Tuple{DID{Int}}) == false
    @test isparameterbound(Coupling, :iids, ID{DID{Int}}) == true
    @test isparameterbound(Coupling, :subscripts, Subscripts{Tuple{NamedTuple{(:nambu,), Tuple{Subscript{Tuple{Int}}}}}}) == true
    @test isparameterbound(Coupling, :subscripts, Subscripts) == true
    @test isparameterbound(Coupling, :subscripts, Subscripts{Tuple{}, Tuple{}}) == false

    tc = DCoupling(2.0, (2,))
    @test id(tc) == ID(CompositeIID(tc.iids), tc.subscripts)
    @test rank(tc) == rank(typeof(tc)) == 1
    @test tc == Coupling(2.0, id(tc))
    @test ID{SimpleIID}(tc) == tc.iids
    @test Subscripts(tc) == tc.subscripts

    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>DFock(2))
    @test couplingcenters(tc, bond, Val(:Mu)) == (1,)
    @test couplingpoints(tc, bond, Val(:Mu)) == (point,)
    @test couplinginternals(tc, bond, hilbert, Val(:Mu)) == (DFock(2),)

    tc₁ = DCoupling(1.5, (1, 2))
    tc₂ = DCoupling(2.0, subscript"[a b](a < b)")
    @test tc₁*tc₂ == Coupling(3.0, ID(DID(1), DID(2), DID(:a), DID(:b)), Subscripts((nambu=Subscript((1, 2)),), (nambu=subscript"[a b](a < b)",)))

    ex = expand(tc₁, bond, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CompositeIndex{Index{DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(1.5, ID(
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        ]

    ex = expand(tc₂, bond, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CompositeIndex{Index{DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(2.0, ID(
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
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
    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.site=>DFock(2))
    term = Term{:Mu}(:μ, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3, ishermitian=true)
    @test term|>kind == term|>typeof|>kind == :Mu
    @test term|>id == term|>typeof|>id == :μ
    @test term|>valtype == term|>typeof|>valtype == Float
    @test term|>rank == term|>typeof|>rank == 2
    @test term|>ismodulatable == term|>typeof|>ismodulatable == false
    @test term == deepcopy(term)
    @test isequal(term, deepcopy(term))
    @test repr(term, Bond(point), hilbert) == "Mu: 4.5*ph[2 1]"

    p₁ = Point(1, (0.0, 0.0), (0.0, 0.0))
    p₂ = Point(2, (1.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    term = Term{:Mu}(:μ, 1.5, 0,
        couplings=bond->bond[1].site%2==0 ? Couplings(DCoupling(1.0, (2, 2))) : Couplings(DCoupling(1.0, (1, 1))),
        ishermitian=true,
        amplitude=bond->3,
        modulate=true
    )
    @test term|>ismodulatable == term|>typeof|>ismodulatable == true
    @test repr(term, Bond(p₁), hilbert) == "Mu: 4.5*ph[1 1]"
    @test repr(term, Bond(p₂), hilbert) == "Mu: 4.5*ph[2 2]"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test term.modulate(μ=4.0) == 4.0
    @test isnothing(term.modulate(t=1.0))
    @test update!(term, μ=4.25) == replace(term, value=4.25)
    @test term.value == 4.25

    term = Term{:Hp}(:t, 1.5, 1, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, ishermitian=false)
    @test repr(term, Bond(1, p₁, p₂), hilbert) == "Hp: 4.5*ph[2 1] + h.c."
end

@testset "expand" begin
    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>DFock(2))
    term = Term{:Mu}(:μ, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, ishermitian=true)
    operators = Operators(
        Operator(+2.25,
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    bond = Bond(1, Point(2, (1.5, 1.5), (1.0, 1.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    term = Term{:Hp}(:t, 1.5, 1, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, ishermitian=false)
    operators = Operators(
        Operator(4.5,
            CompositeIndex(Index(2, DID(2)), SVector(1.5, 1.5), SVector(1.0, 1.0)),
            CompositeIndex(Index(1, DID(1)), SVector(0.5, 0.5), SVector(0.0, 0.0))
            )
        )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>DFock(2) for site=1:length(lattice))
    term = Term{:Mu}(:μ, 1.5, 0, couplings=@couplings(DCoupling(1.0, (2, 1))), amplitude=bond->3.0, ishermitian=true)
    operators = Operators(
        Operator(+2.25,
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    @test expand(term, bs, hilbert, half=true) == operators
    @test expand(term, bs, hilbert, half=false) == operators*2
end

@testset "script" begin
    index = CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(1.0, 0.0))
    latex = LaTeX{(:nambu,), (:site,)}('c', vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0)))
    @test script(Val(:rcoordinate), index) == "[0.0, 0.0]"
    @test script(Val(:icoordinate), index) == "[1.0, 0.0]"
    @test script(Val(:integercoordinate), index; vectors=get(latex.options, :vectors, nothing)) == "[1, 0]"
    @test script(Val(:site), index) == 1
    @test script(Val(:nambu), index) == "\\dagger"
end

@testset "Boundary" begin
    op = Operator(4.5,
        CompositeIndex(Index(1, DID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
        CompositeIndex(Index(2, DID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
        )
    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    M = reparameter(typeof(op), :value, Complex{Float64})
    @test valtype(typeof(bound), typeof(op)) == M
    @test keys(bound) == keys(typeof(bound)) == (:θ₁, :θ₂)
    @test bound == deepcopy(bound)
    @test bound≠Boundary{(:ϕ₁, :ϕ₂)}(bound.values, bound.vectors)
    @test isequal(bound, deepcopy(bound))
    @test !isequal(bound, Boundary{(:ϕ₁, :ϕ₂)}(bound.values, bound.vectors))

    another = Boundary{(:θ₁, :θ₂)}([0.0, 0.0], [[2.0, 0.0], [0.0, 2.0]])
    @test merge!(deepcopy(bound), another) == another
    @test replace(bound; values=another.values, vectors=another.vectors) == another

    @test bound(op) ≈ replace(op, 4.5*exp(2im*pi*0.3))
    @test bound(op, origin=[0.05, 0.15]) ≈ replace(op, 4.5*exp(2im*pi*0.1))
    update!(bound, θ₁=0.3)
    @test bound(op) ≈ replace(op, 4.5*exp(2im*pi*0.5))
    @test bound(op, origin=[0.1, 0.1]) ≈ replace(op, 4.5*exp(2im*pi*0.3))

    ops = Operators{M}(op)
    @test valtype(typeof(bound), typeof(ops)) == typeof(ops)
    @test bound(ops) ≈ Operators(replace(op, 4.5*exp(2im*pi*0.5)))
    @test map!(bound, ops) ≈ ops ≈ Operators(replace(op, 4.5*exp(2im*pi*0.5)))

    @test valtype(typeof(plain), typeof(op)) == typeof(op)
    @test valtype(typeof(plain), typeof(ops)) == typeof(ops)
    @test plain(op) == op
    @test plain(ops) == ops
    @test update!(plain) == plain
    @test replace(plain; values=another.values, vectors=another.vectors) == plain
end
