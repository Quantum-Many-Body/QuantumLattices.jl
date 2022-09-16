using LaTeXStrings: latexstring
using LinearAlgebra: dot
using Printf: @printf, @sprintf
using QuantumLattices.Essentials: kind, reset!, update!
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.QuantumOperators: ID, LaTeX, Operator, Operators, id, idtype, latexformat, sequence
using QuantumLattices.Essentials.Spatials: Bond, Lattice, Point, bonds, decompose, icoordinate, rcoordinate
using QuantumLattices.Interfaces: ⊕, ⊗, dimension, expand, rank
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Traits: contentnames, getcontent, isparameterbound, parameternames, reparameter
using StaticArrays: SVector

import QuantumLattices.Essentials.DegreesOfFreedom: Constraint, iidtype, isconcreteiid, statistics
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
@inline isconcreteiid(::Type{DID{Int}}) = true
@inline iidtype(::Type{DID}, ::Type{T}) where {T<:Union{Int, Symbol}} = DID{T}
@inline Constraint(iids::NTuple{N, DID{Int}}) where N = Constraint{N}()

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

@inline script(::Val{:nambu}, did::DID; kwargs...) = did.nambu==2 ? "\\dagger" : did.nambu==1 ? "" : did.nambu

@inline latexname(::Type{<:CompositeIndex{<:Index{<:DID}}}) = Symbol("CompositeIndex{Index{DID}}")
latexformat(CompositeIndex{<:Index{<:DID}}, LaTeX{(:nambu,), (:site,)}('d'))
latexformat(DID, LaTeX{(:nambu,), ()}('d'))

const DCoupling{V, I<:ID{DID}, C<:Constraint} = Coupling{V, I, C}
@inline DCoupling(value, nambus::NTuple{N, Int}) where N = Coupling(value, ID(DID, nambus))

@testset "IID" begin
    did = DID(1)
    @test statistics(did) == statistics(typeof(did)) == :f
    @test isconcreteiid(did) == isconcreteiid(typeof(did)) == true
    @test isconcreteiid(DID(:a)) == false

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
    @test rcoordinate(opt) == SVector(-0.5, 0.0)
    @test icoordinate(opt) == SVector(-1.0, 0.0)

    opt = Operator(1.0, CompositeIndex(Index(1, DID(2)), SVector(0.5, 0.0), SVector(1.0, 0.0)))
    @test rcoordinate(opt) == SVector(0.5, 0.0)
    @test icoordinate(opt) == SVector(1.0, 0.0)
end

@testset "Hilbert" begin
    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    @test hilbert[1] == hilbert[2] == DFock(2)

    hilbert = Hilbert(1=>DFock(2), 2=>DFock(3))
    @test hilbert[1]==DFock(2) && hilbert[2]==DFock(3)
end

@testset "IIDSpace" begin
    DID₁, DID₂, it = DID(2), DID(:σ), DFock(2)
    iidspace = IIDSpace(DID₁⊗DID₂, it⊗it)
    @test eltype(iidspace) == eltype(typeof(iidspace)) == CompositeIID{Tuple{DID{Int}, DID{Int}}}
    @test length(iidspace) == 2
    @test iidspace[1]==DID(2)⊗DID(1) && iidspace[2]==DID(2)⊗DID(2)
    @test expand((DID₁, DID₂), (it, it)) == iidspace
    @test collect(iidspace) == [DID(2)⊗DID(1), DID(2)⊗DID(2)]
end

@testset "Constraint" begin
    constraint = Constraint(DID(:a), DID(:a))
    @test deepcopy(constraint) == constraint
    @test isequal(deepcopy(constraint), constraint)
    @test hash(constraint) == hash(((2,), constraint.representations))
    @test rank(constraint) == rank(typeof(constraint)) == 2
    @test rank(constraint, 1) == rank(typeof(constraint), 1) == 2
    @test match(constraint, (DID(2), DID(2)))
    @test !match(constraint, (DID(2), DID(1)))

    constraint = (@iids DID(a) DID(a) DID(b))[2]
    @test match(constraint, (DID(1), DID(1), DID(2)))
    @test !match(constraint, (DID(1), DID(2), DID(2)))

    another = @iids(DID(2), DID(1))[2]
    @test match(another, (DID(2), DID(1)))
    @test match(another, (DID(2), DID(2)))

    another = @iids(DID(a), DID(b); constraint=a<b)[2]
    @test match(another, (DID(1), DID(2)))
    @test !match(another, (DID(2), DID(1)))

    product = constraint*another
    @test match(product, (DID(2), DID(2), DID(3), DID(1), DID(2)))
    @test !match(product, (DID(1), DID(2), DID(3), DID(1), DID(2)))
    @test !match(product, (DID(2), DID(2), DID(3), DID(2), DID(1)))
    @test !match(product, (DID(2), DID(1), DID(3), DID(2), DID(1)))

    constraint = Constraint{2}(Diagonal(:nambu))
    @test match(constraint, (DID(2), DID(2)))
    @test !match(constraint, (DID(2), DID(1)))
end

@testset "Coupling" begin
    @test parameternames(Coupling) == (:value, :iids, :constraint)
    @test isparameterbound(Coupling, :iids, Tuple{DID{Int}}) == false
    @test isparameterbound(Coupling, :iids, ID{DID{Int}}) == true
    @test isparameterbound(Coupling, :constraint, Constraint) == true
    @test isparameterbound(Coupling, :constraint, Constraint{(2, 2)}) == true
    @test isparameterbound(Coupling, :constraint, Constraint{(2,), 1, Tuple{Pattern}}) == false

    tc = DCoupling(2.0, (2,))
    @test id(tc) == (tc.iids, tc.constraint)
    @test idtype(typeof(tc)) == typeof(id(tc))
    @test rank(tc) == rank(typeof(tc)) == 1
    @test Coupling(2.0, id(tc)) == tc
    @test CompositeIID(tc) == CompositeIID(tc.iids)
    @test Constraint(tc) == tc.constraint
    @test length(tc) == length(typeof(tc)) == 1
    @test eltype(tc) == eltype(typeof(tc)) == typeof(tc)
    @test collect(tc) == [tc]

    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>DFock(2))
    @test sitestructure(Val(:Mu), Val(rank(tc)), length(bond)) == (1,)

    tc₁ = DCoupling(1.5, (1, 2))
    tc₂ = Coupling(2.0, @iids(DID(a), DID(b); constraint=a<b))
    tc =  tc₁*tc₂
    @test tc₁*tc₂ == Coupling(
        3.0,
        (DID(1), DID(2), DID(:a), DID(:b)),
        Constraint{(2, 2)}((tc₁.constraint.representations[1], tc₂.constraint.representations[1]), (tc₁.constraint.conditions[1], tc₂.constraint.conditions[1]))
    )
    @test string(tc)=="3.0 [DID(1) DID(2)] ⋅ ∑[DID(a) DID(b)](a < b)"
    @test latexstring(tc) == "3.0 d^{}_{} d^{\\dagger}_{} \\cdot \\sum_{a < b} d^{a}_{} d^{b}_{}"

    ex = expand(Val(:term), tc₁, bond, hilbert)
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CompositeIndex{Index{DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(1.5, ID(
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        ]

    ex = expand(Val(:term), tc₂, bond, hilbert)
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CompositeIndex{Index{DID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(2.0, ID(
            CompositeIndex(Index(1, DID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CompositeIndex(Index(1, DID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        ]
end

@testset "MatrixCoupling" begin
    component = Component([1, 2], [2, 1], [-1 0; 0 1])
    @test dimension(component) == 2
    @test component[1] == (1, 2, -1)
    @test component[2] == (2, 1, +1)

    mc = MatrixCoupling{DID}(component)
    @test constrainttype(typeof(mc)) == Constraint{(2,), 1, Tuple{Pattern}}
    @test eltype(typeof(mc)) == Coupling{Int, Tuple{DID{Int}, DID{Int}}, Constraint{(2,), 1, Tuple{Pattern}}}
    @test mc[1] == Coupling(-1, DID(1), DID(2))
    @test mc[2] == Coupling(+1, DID(2), DID(1))

    another = MatrixCoupling{DID}(Component([:α], [:α], hcat(2.0)))
    @test constrainttype(typeof(another)) == Constraint{(2,), 1, Tuple{Pattern}}
    @test another[1] == Coupling(2.0, DID(:α), DID(:α))

    mcs = mc*another
    @test mcs == MatrixCouplingProd(mc, another)
    @test eltype(mcs) == Coupling{Float64, Tuple{DID{Int}, DID{Int}, DID{Symbol}, DID{Symbol}}, Constraint{(2, 2), 2, Tuple{Pattern, Pattern}}}
    @test mcs[1] == Coupling(-2.0, (DID(1), DID(2), DID(:α), DID(:α)), Constraint{(2, 2)}(("pattern", "pattern"), (noconstrain, Pattern(DID(:α), DID(:α)))))
    @test mcs[2] == Coupling(+2.0, (DID(2), DID(1), DID(:α), DID(:α)), Constraint{(2, 2)}(("pattern", "pattern"), (noconstrain, Pattern(DID(:α), DID(:α)))))

    @test mcs*mc == MatrixCouplingProd(mc, another, mc)
    @test mc*mcs == MatrixCouplingProd(mc, mc, another)
    @test mcs*mcs == MatrixCouplingProd(mc, another, mc, another)

    mc₁ = MatrixCoupling{DID}(Component([1, 2], [2, 1], [-1 0; 0 1]))
    mc₂ = MatrixCoupling{DID}(Component([1, 2], [2, 1], [0 1; 1 0]))
    mcs = mc₁ + mc₂
    @test mcs == MatrixCouplingSum(mc₁, mc₂)
    @test eltype(mcs) == Coupling{Int, Tuple{DID{Int}, DID{Int}}, Constraint{(2,), 1, Tuple{Pattern}}}
    @test collect(mcs) == [
        Coupling(-1, DID(1), DID(2)), Coupling(+1, DID(2), DID(1)),
        Coupling(+1, DID(2), DID(2)), Coupling(+1, DID(1), DID(1))
    ]

    @test mcs+mc₁ == MatrixCouplingSum(mc₁, mc₂, mc₁)
    @test mc₁+mcs == MatrixCouplingSum(mc₁, mc₁, mc₂)
    @test mcs+mcs == MatrixCouplingSum(mc₁, mc₂, mc₁, mc₂)
end

@testset "TermFunction" begin
    bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))

    ta = TermAmplitude()
    @test ta(bond) == 1
    ta = TermAmplitude(bond::Bond->bond.kind+3)
    @test ta(bond) == 4

    tcs = DCoupling(1.0, (1, 1)) + DCoupling(2.0, (2, 2))
    termcouplings = TermCoupling(tcs)
    @test termcouplings == deepcopy(TermCoupling(tcs))
    @test isequal(termcouplings, deepcopy(TermCoupling(tcs)))
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == eltype(typeof(tcs))
    @test termcouplings(bond) == tcs

    bond₁ = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
    bond₂ = Bond(2, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))

    fx = bond::Bond -> bond.kind==1 ? DCoupling(1.0, (1, 1)) : DCoupling(1.0, (2, 2))
    termcouplings = TermCoupling(fx)
    @test termcouplings == TermCoupling{typejoin(typeof(fx(bond₁)), typeof(fx(bond₂)))}(fx)
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typejoin(typeof(fx(bond₁)), typeof(fx(bond₂)))
    @test termcouplings(bond₁) == fx(bond₁)
    @test termcouplings(bond₂) == fx(bond₂)

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
    term = Term{:Mu}(:μ, 1.5, 0, DCoupling(1.0, (2, 1)), true; amplitude=bond->3, modulate=false)
    @test term|>kind == term|>typeof|>kind == :Mu
    @test term|>id == term|>typeof|>id == :μ
    @test term|>valtype == term|>typeof|>valtype == Float
    @test term|>rank == term|>typeof|>rank == 2
    @test term|>ismodulatable == term|>typeof|>ismodulatable == false
    @test term == deepcopy(term)
    @test isequal(term, deepcopy(term))
    @test repr(term, Bond(point), hilbert) == "Mu: 4.5 DID(2) DID(1)"

    p₁ = Point(1, (0.0, 0.0), (0.0, 0.0))
    p₂ = Point(2, (1.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    term = Term{:Mu}(:μ, 1.5, 0, bond->bond[1].site%2==0 ? DCoupling(1.0, (2, 2)) : DCoupling(1.0, (1, 1)), true; amplitude=bond->3)
    @test term|>ismodulatable == term|>typeof|>ismodulatable == true
    @test repr(term, Bond(p₁), hilbert) == "Mu: 4.5 DID(1) DID(1)"
    @test repr(term, Bond(p₂), hilbert) == "Mu: 4.5 DID(2) DID(2)"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test term.modulate(μ=4.0) == 4.0
    @test isnothing(term.modulate(t=1.0))
    @test update!(term, μ=4.25) == replace(term, value=4.25)
    @test term.value == 4.25

    term = Term{:Hp}(:t, 1.5, 1, DCoupling(1.0, (2, 1)), false; amplitude=bond->3.0)
    @test repr(term, Bond(1, p₁, p₂), hilbert) == "Hp: 4.5 DID(2) DID(1) + h.c."
end

@testset "expand" begin
    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>DFock(2))
    term = Term{:Mu}(:μ, 1.5, 0, DCoupling(1.0, (2, 1)), true; amplitude=bond->3.0)
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
    term = Term{:Hp}(:t, 1.5, 1, DCoupling(1.0, (2, 1)), false; amplitude=bond->3.0)
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
    term = Term{:Mu}(:μ, 1.5, 0, DCoupling(1.0, (2, 1)), true; amplitude=bond->3.0)
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
    @test script(Val(:site), index) == "1"
    @test script(Val(:nambu), index) == "\\dagger"
end

@testset "Metric" begin
    valtype(Metric, AbstractOID) = AbstractOID

    m = OperatorUnitToTuple((:site, :nambu))
    @test m == OperatorUnitToTuple(:site, :nambu)
    @test isequal(m, OperatorUnitToTuple(:site, :nambu))
    @test keys(m) == keys(typeof(m)) == (:site, :nambu)
    @test OperatorUnitToTuple(Index{DID{Int}}) == OperatorUnitToTuple(:site, :nambu)
    @test OperatorUnitToTuple(Hilbert{DFock}) == OperatorUnitToTuple(:site, :nambu)

    index = CompositeIndex(Index(4, DID(1)), SVector(0.5, 0.0), SVector(1.0, 0.0))
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
