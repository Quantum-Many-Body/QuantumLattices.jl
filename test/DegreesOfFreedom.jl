using LaTeXStrings: latexstring
using LinearAlgebra: dot, ishermitian
using Printf: @printf
using QuantumLattices: ⊕, ⊗, expand, kind, rank, reset!, update!
using QuantumLattices.DegreesOfFreedom
using QuantumLattices.QuantumOperators: ID, LaTeX, Operator, Operators, id, latexformat, sequence
using QuantumLattices.Spatials: Bond, Point, decompose, icoordinate, rcoordinate
using QuantumLattices.Toolkit: Float, contentnames, parameternames, reparameter
using StaticArrays: SVector

import QuantumLattices: permute
import QuantumLattices.DegreesOfFreedom: indextype, isdefinite, statistics
import QuantumLattices.QuantumOperators: latexname, script
import QuantumLattices.Toolkit: shape

struct DID{N<:Union{Int, Symbol, Colon}} <: SimpleInternalIndex
    nambu::N
end
@inline Base.show(io::IO, did::DID) = @printf io "DID(%s)" did.nambu==Colon() ? ":" : string(did.nambu)
@inline Base.show(io::IO, ::Type{<:DID}) = @printf io "%s" "DID"
@inline Base.adjoint(sl::DID{Int}) = DID(3-sl.nambu)
@inline statistics(::Type{<:DID}) = :f
function permute(did₁::DID, did₂::DID)
    @assert did₁ ≠ did₂ "permute error: two identical fermionic indexes should vanish due to the fermionic statistics."
    return (Operator(1), Operator(-1, did₂, did₁))
end
@inline isdefinite(::Type{DID{Int}}) = true
@inline script(did::DID, ::Val{:nambu}; kwargs...) = did.nambu==Colon() ? ":" : string(did.nambu)
function Base.angle(id::CoordinatedIndex{Index{DID{Int}, Int}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase = (length(vectors) == 1) ? 2pi*dot(decompose(id.icoordinate, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return (id.index.internal.nambu == 1) ? phase : -phase
end
@inline indextype(::Type{DID}, ::Type{T}) where {T<:Union{Int, Symbol, Colon}} = DID{T}

struct DFock <: SimpleInternal{DID{Int}}
    nnambu::Int
end
@inline shape(f::DFock) = (1:f.nnambu,)
@inline Base.show(io::IO, ::Type{DFock}) = @printf io "%s" "DFock"
@inline Base.convert(::Type{<:DID}, i::CartesianIndex, ::DFock) = DID(i.I...)
@inline Base.convert(::Type{<:CartesianIndex}, did::DID{Int}, ::DFock) = CartesianIndex(did.nambu)
@inline shape(internal::DFock, ::DID) = (1:internal.nnambu,)
@inline shape(::DFock, index::DID{Int}) = (index.nambu:index.nambu,)
@inline latexname(::Type{<:CoordinatedIndex{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}}) = Symbol("CoordinatedIndex{Index{DID, Union{Int, Ordinal, Colon}}}")
@inline latexname(::Type{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}) = Symbol("Index{DID, Union{Int, Ordinal, Colon}}")
latexformat(CoordinatedIndex{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}, LaTeX{(), (:site, :nambu)}('d'))
latexformat(Index{<:DID, <:Union{Int, Ordinal, Colon}}, LaTeX{(), (:site, :nambu)}('d'))
latexformat(DID, LaTeX{(), (:nambu,)}('d'))

@testset "InternalIndex" begin
    did = DID(1)
    @test statistics(did) == statistics(typeof(did)) == :f
    @test isdefinite(did) == isdefinite(typeof(did)) == true
    @test isdefinite(DID(:a)) == isdefinite(DID{Symbol}) == false

    did₁, did₂ = DID(1), DID(2)
    ciid = InternalIndexProd(did₁, did₂)
    @test string(ciid) == "DID(1) ⊗ DID(2)"
    @test length(ciid) == rank(ciid) == rank(typeof(ciid)) == 2
    @test ciid[1]==ciid[begin]==did₁ && ciid[2]==ciid[end]==did₂
    @test ciid[2:2] == InternalIndexProd(did₂)
    @test ciid.contents==(DID(1), DID(2)) && ciid.nambus==(1, 2)
    @test isdefinite(ciid)
    @test did₁⊗did₂ == ciid
    @test did₁⊗ciid == InternalIndexProd(did₁, did₁, did₂)
    @test ciid⊗did₁ == InternalIndexProd(did₁, did₂, did₁)
    @test ciid⊗ciid == InternalIndexProd(did₁, did₂, did₁, did₂)
end

@testset "SimpleInternal" begin
    it = DFock(2)
    @test string(it) == "DFock(nnambu=2)"
    @test collect(it) == [DID(1), DID(2)]
    @test statistics(it) == statistics(typeof(it)) == :f
    @test match(DID(1), it) && match(DID, it) && match(DID(1), DFock) && match(DID, DFock)
    @test filter(DID(1), it) == filter(DID, it) == it
    @test filter(DID(1), DFock) == filter(DID, DFock) == DFock
end

@testset "CompositeInternal" begin
    it₁, it₂ = DFock(2), DFock(3)

    ci = InternalSum(it₁, it₂)
    @test eltype(ci) == eltype(typeof(ci)) == DID{Int}
    @test rank(ci) == rank(typeof(ci)) == 2
    @test string(ci) == "DFock(nnambu=2) ⊕ DFock(nnambu=3)"
    @test collect(ci) == [it₁[1], it₁[2], it₂[1], it₂[2], it₂[3]]
    @test it₁⊕it₂ == ci
    @test it₁⊕ci == InternalSum(it₁, it₁, it₂)
    @test ci⊕it₁ == InternalSum(it₁, it₂, it₁)
    @test ci⊕ci == InternalSum(it₁, it₂, it₁, it₂)
    @test filter(DID(1), ci) == filter(DID, ci) == ci
    @test filter(DID(1), typeof(ci)) == filter(DID, typeof(ci)) == typeof(ci)

    ci = InternalProd(it₁, it₂)
    @test eltype(ci) == eltype(typeof(ci)) == InternalIndexProd{Tuple{DID{Int}, DID{Int}}}
    @test rank(ci) == rank(typeof(ci)) == 2
    @test string(ci) == "DFock(nnambu=2) ⊗ DFock(nnambu=3)"
    @test collect(ci) == InternalIndexProd.([(it₁[1], it₂[1]), (it₁[2], it₂[1]), (it₁[1], it₂[2]), (it₁[2], it₂[2]), (it₁[1], it₂[3]), (it₁[2], it₂[3])])
    @test it₁⊗it₂ == ci
    @test it₁⊗ci == InternalProd(it₁, it₁, it₂)
    @test ci⊗it₁ == InternalProd(it₁, it₂, it₁)
    @test ci⊗ci == InternalProd(it₁, it₂, it₁, it₂)
    @test filter(DID(1), ci) == filter(DID, ci) == ci
    @test filter(DID(1), typeof(ci)) == filter(DID, typeof(ci)) == typeof(ci)
end

@testset "ConstrainedInternal" begin
    allequal = AllEqual(DID{Int})
    @test string(allequal) == "AllEqual()"
    @test allequal(DID(1) ⊗ DID(1)) && allequal(DID(1) ⊗ DID(2))

    allequal = AllEqual(DID{Colon})
    @test string(allequal) == "AllEqual(:nambu)"
    @test allequal(DID(1) ⊗ DID(1)) && !allequal(DID(1) ⊗ DID(2))

    pattern = InternalPattern(DID(:)⊗DID(:), allequal, "AllEqual(:nambu)")
    @test pattern == InternalPattern(DID(:)⊗DID(:))
    @test isequal(pattern, InternalPattern(DID(:)⊗DID(:)))
    @test hash(pattern, UInt(8)) == hash((2, DID(:), DID(:), "AllEqual(:nambu)"), UInt(8))
    @test string(pattern) == "∑[DID(:) DID(:)]"
    @test partition(pattern) == partition(typeof(pattern)) == (2,)
    @test rank(pattern) == rank(typeof(pattern)) == 2
    @test rank(pattern, 1) == rank(typeof(pattern), 1) == 2
    @test match(pattern, DID(1)⊗DID(1)) && !match(pattern, DID(1)⊗DID(2))
    @test latexstring(pattern) == "\\sum_{} d^{}_{:} d^{}_{:}"

    another = pattern ⊗ pattern
    @test string(another) == "∑[DID(:) DID(:)] ⊗ ∑[DID(:) DID(:)]"
    @test partition(another) == partition(typeof(another)) == (2, 2)
    @test rank(another) == rank(typeof(another)) == 4
    @test rank(another, 1) == rank(typeof(another), 1) == 2
    @test rank(another, 2) == rank(typeof(another), 2) == 2
    @test match(another, DID(1)⊗DID(1)⊗DID(2)⊗DID(2)) && !match(another, DID(1)⊗DID(1)⊗DID(2)⊗DID(1))
    @test latexstring(another) == "\\sum_{} d^{}_{:} d^{}_{:} \\cdot \\sum_{} d^{}_{:} d^{}_{:}"

    pattern = InternalPattern(DID(:a)⊗DID(:b), index->index[1].nambu<index[2].nambu, "a < b")
    @test string(pattern) == "∑[DID(a) DID(b)](a < b)"
    @test match(pattern, DID(1)⊗DID(2)) && !match(pattern, DID(2)⊗DID(1))
    @test latexstring(pattern) == "\\sum_{a < b} d^{}_{a} d^{}_{b}"

    another = pattern ⊗ pattern
    @test string(another) == "∑[DID(a) DID(b)](a < b) ⊗ ∑[DID(a) DID(b)](a < b)"
    @test match(another, DID(1)⊗DID(2)⊗DID(3)⊗DID(4)) && !match(another, DID(2)⊗DID(1)⊗DID(3)⊗DID(4)) && !match(another, DID(1)⊗DID(2)⊗DID(4)⊗DID(3))
    @test latexstring(another) == "\\sum_{a < b} d^{}_{a} d^{}_{b} \\cdot \\sum_{a < b} d^{}_{a} d^{}_{b}"

    con = ConstrainedInternal(DFock(2), InternalPattern(DID(:)))
    @test eltype(con) == eltype(typeof(con)) == InternalIndexProd{Tuple{DID{Int}}}
    @test collect(con) == InternalIndexProd.([DID(1), DID(2)])

    con = ConstrainedInternal(DFock(2)⊗DFock(2), InternalPattern(DID(:)⊗DID(:)))
    @test eltype(con) == eltype(typeof(con)) == InternalIndexProd{Tuple{DID{Int}, DID{Int}}}
    @test collect(con) == [DID(1)⊗DID(1), DID(2)⊗DID(2)]

    con = ConstrainedInternal(DFock(2)⊗DFock(2)⊗DFock(2)⊗DFock(2), InternalPattern(DID(:)⊗DID(:))⊗InternalPattern(DID(:)⊗DID(:)))
    @test eltype(con) == eltype(typeof(con)) == InternalIndexProd{Tuple{DID{Int}, DID{Int}, DID{Int}, DID{Int}}}
    @test collect(con) == [DID(1)⊗DID(1)⊗DID(1)⊗DID(1), DID(2)⊗DID(2)⊗DID(1)⊗DID(1), DID(1)⊗DID(1)⊗DID(2)⊗DID(2), DID(2)⊗DID(2)⊗DID(2)⊗DID(2)]
end

@testset "Index" begin
    @test (4, 3, 2, 1)[1ˢᵗ] == 4
    @test (4, 3, 2, 1)[2ⁿᵈ] == 3
    @test (4, 3, 2, 1)[3ʳᵈ] == 2
    @test (4, 3, 2, 1)[4ᵗʰ] == 1

    @test parameternames(Index) == (:internal, :site)

    index = Index(4, DID(1))
    @test indextype(index) == indextype(typeof(index)) == DID{Int}
    @test index' == Index{DID}(4, 2)
    @test statistics(index) == statistics(typeof(index)) == :f
    @test ishermitian(ID(index', index)) == true
    @test ishermitian(ID(index, index)) == false
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test isdefinite((index, index)) == isdefinite(typeof((index, index))) == true

    @test string(Index{DID}(:, 2)) == "Index{DID}(:, 2)"
    @test string(Index{DID}(1ˢᵗ, 2)) == "Index{DID}(1ˢᵗ, 2)"
    @test string(Index{DID}(2ⁿᵈ, 2)) == "Index{DID}(2ⁿᵈ, 2)"
    @test string(Index{DID}(3ʳᵈ, 2)) == "Index{DID}(3ʳᵈ, 2)"
    @test string(Index{DID}(4ᵗʰ, 2)) == "Index{DID}(4ᵗʰ, 2)"

    @test script(Index{DID}(1, 2), Val(:site)) == "1"
    @test script(Index{DID}(:, 2), Val(:site)) == ":"
    @test script(Index{DID}(1, 2), Val(:nambu)) == "2"

    index₁, index₂ = Index{DID}(1, 2), Index{DID}(1, 1)
    @test permute(index₁, index₂) == (Operator(1), Operator(-1, index₂, index₁))

    index₁, index₂ = Index{DID}(1, 2), Index{DID}(2, 2)
    @test permute(index₁, index₂) == (Operator(-1, index₂, index₁),)

    @test indextype(DFock) == Index{DID{Int}, Int}
end

@testset "CoordinatedIndex" begin
    @test contentnames(CompositeIndex) == (:index,)
    @test parameternames(CompositeIndex) == (:index,)

    @test contentnames(CoordinatedIndex) == (:index, :rcoordinate, :icoordinate)
    @test parameternames(CoordinatedIndex) == (:index, :coordination)

    index = CoordinatedIndex(Index{DID}(1, 1), [0.0, -0.0], [0.0, 0.0])
    @test indextype(index) == indextype(typeof(index)) == Index{DID{Int}, Int}
    @test statistics(index) == statistics(typeof(index)) == :f
    @test hash(index, UInt(1)) == hash(CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.0), SVector(0.0, 1.0)), UInt(1))
    @test propertynames(ID(index)) == (:indexes, :rcoordinates, :icoordinates)
    @test string(index) == "CoordinatedIndex(Index{DID}(1, 1), [0.0, 0.0], [0.0, 0.0])"
    @test index' == CoordinatedIndex(Index{DID}(1, 2), rcoordinate=SVector(0.0, 0.0), icoordinate=SVector(0.0, 0.0))
    @test ID(index', index)' == ID(index', index)
    @test ishermitian(ID(index', index)) && !ishermitian(ID(index, index))

    index = CoordinatedIndex(Index{DID}(1, 2), SVector(0.0, 0.0), SVector(1.0, 0.0))
    @test script(index, Val(:rcoordinate)) == "[0.0, 0.0]"
    @test script(index, Val(:icoordinate)) == "[1.0, 0.0]"
    @test script(index, Val(:integercoordinate); vectors=(SVector(1.0, 0.0), SVector(0.0, 1.0))) == "[1, 0]"
    @test script(index, Val(:site)) == "1"
    @test script(index, Val(:nambu)) == "2"

    index₁ = CoordinatedIndex(Index{DID}(1, 1), (1.0, 0.0), (0.0, 0.0))
    index₂ = CoordinatedIndex(Index{DID}(1, 1), (2.0, 0.0), (0.0, 0.0))
    @test permute(index₁, index₂) == (Operator(-1, index₂, index₁),)

    index₁ = CoordinatedIndex(Index{DID}(1, 1), (1.0, 0.0), (0.0, 0.0))
    index₂ = CoordinatedIndex(Index{DID}(1, 2), (1.0, 0.0), (0.0, 0.0))
    @test permute(index₁, index₂) == (Operator(1), Operator(-1, index₂, index₁),)

    @test indextype(DFock, Point{2, Float}) == CoordinatedIndex{Index{DID{Int}, Int}, SVector{2, Float}}

    op = Operator(1.0,
        CoordinatedIndex(Index{DID}(1, 2), SVector(0.5, 0.5), SVector(1.0, 1.0)),
        CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.5), SVector(0.0, 1.0))
        )
    @test rcoordinate(op) == SVector(-0.5, 0.0)
    @test icoordinate(op) == SVector(-1.0, 0.0)

    op = Operator(1.0, CoordinatedIndex(Index{DID}(1, 2), SVector(0.5, 0.0), SVector(1.0, 0.0)))
    @test rcoordinate(op) == SVector(0.5, 0.0)
    @test icoordinate(op) == SVector(1.0, 0.0)
end

@testset "Hilbert" begin
    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    @test hilbert[1] == hilbert[2] == DFock(2)
    @test hilbert == Hilbert(DFock(2), 2)

    hilbert = Hilbert(1=>DFock(2), 2=>DFock(3))
    @test hilbert[1]==DFock(2) && hilbert[2]==DFock(3)
    @test hilbert == Hilbert(DFock(2), DFock(3)) == Hilbert([DFock(2), DFock(3)])
end

@testset "Pattern" begin
    @test parameternames(Pattern) == (:sites, :internal)

    pattern = @pattern(Index(1ˢᵗ, DID(a)), Index(1ˢᵗ, DID(a)), Index(2ⁿᵈ, DID(b)), Index(2ⁿᵈ, DID(b)))
    @test pattern == @pattern(Index{DID}(1ˢᵗ, a), Index{DID}(1ˢᵗ, a), Index{DID}(2ⁿᵈ, b), Index{DID}(2ⁿᵈ, b))
    @test isequal(pattern, @pattern(Index{DID}(1ˢᵗ, a), Index{DID}(1ˢᵗ, a), Index{DID}(2ⁿᵈ, b), Index{DID}(2ⁿᵈ, b)))
    @test hash(pattern, UInt(8)) == hash((pattern.sites, pattern.internal), UInt(8))
    @test string(pattern) == "∑[Index{DID}(1ˢᵗ, a) Index{DID}(1ˢᵗ, a) Index{DID}(2ⁿᵈ, b) Index{DID}(2ⁿᵈ, b)]"
    @test rank(pattern) == rank(typeof(pattern)) == 4
    @test latexstring(pattern) == "\\sum_{} d^{}_{1ˢᵗ,\\,a} d^{}_{1ˢᵗ,\\,a} d^{}_{2ⁿᵈ,\\,b} d^{}_{2ⁿᵈ,\\,b}"
    @test match(pattern.internal, DID(3)⊗DID(3)⊗DID(1)⊗DID(1))
    @test !match(pattern.internal, DID(3)⊗DID(3)⊗DID(1)⊗DID(2))
    @test !match(pattern.internal, DID(3)⊗DID(1)⊗DID(2)⊗DID(2))

    pattern = @pattern(Index{DID}(1ˢᵗ, 1), Index{DID}(1ˢᵗ, a), Index{DID}(2ⁿᵈ, 2), Index{DID}(2ⁿᵈ, b); constraint=a<b)
    @test string(pattern) == "∑[Index{DID}(1ˢᵗ, 1) Index{DID}(1ˢᵗ, a) Index{DID}(2ⁿᵈ, 2) Index{DID}(2ⁿᵈ, b)](a < b)"
    @test latexstring(pattern) == "\\sum_{a < b} d^{}_{1ˢᵗ,\\,1} d^{}_{1ˢᵗ,\\,a} d^{}_{2ⁿᵈ,\\,2} d^{}_{2ⁿᵈ,\\,b}"
    @test match(pattern.internal, DID(1)⊗DID(3)⊗DID(2)⊗DID(4))
    @test !match(pattern.internal, DID(10)⊗DID(3)⊗DID(2)⊗DID(4))
    @test !match(pattern.internal, DID(1)⊗DID(3)⊗DID(2)⊗DID(3))
    @test !match(pattern.internal, DID(1)⊗DID(3)⊗DID(4)⊗DID(5))

    pattern = Pattern(Index{DID}(1ˢᵗ, :), Index{DID}(1ˢᵗ, :))
    @test string(pattern) == "∑[Index{DID}(1ˢᵗ, :) Index{DID}(1ˢᵗ, :)]"
    @test latexstring(pattern) == "\\sum_{} d^{}_{1ˢᵗ,\\,:} d^{}_{1ˢᵗ,\\,:}"
    @test match(pattern.internal, DID(1)⊗DID(1))
    @test !match(pattern.internal, DID(1)⊗DID(3))

    another = pattern ⊗ pattern
    @test string(another) == "∑[Index{DID}(1ˢᵗ, :) Index{DID}(1ˢᵗ, :)] ⊗ ∑[Index{DID}(1ˢᵗ, :) Index{DID}(1ˢᵗ, :)]"
    @test latexstring(another) == "\\sum_{} d^{}_{1ˢᵗ,\\,:} d^{}_{1ˢᵗ,\\,:} \\cdot \\sum_{} d^{}_{1ˢᵗ,\\,:} d^{}_{1ˢᵗ,\\,:}"
    @test match(another.internal, DID(1)⊗DID(1)⊗DID(2)⊗DID(2))
    @test !match(another.internal, DID(1)⊗DID(3)⊗DID(2)⊗DID(2))
    @test !match(another.internal, DID(1)⊗DID(1)⊗DID(3)⊗DID(2))
end

@testset "patternrule" begin
    @test patternrule((1, 2, 3, 4), Val(:)) == (1, 2, 3, 4)
    @test patternrule(DID(1)⊗DID(:), Val(:)) == DID(1)⊗DID(:)
    @test patternrule(DID(1)⊗DID(2), Val(:)) == DID(1)⊗DID(2)
    @test patternrule((1, 2), Val(:), DID, Val(:nambu)) == (1, 2)
    @test patternrule((:, :), Val(:), 1) == (1ˢᵗ, 1ˢᵗ)
    @test patternrule((:, :), Val(:), 2) == (1ˢᵗ, 2ⁿᵈ)
    @test patternrule((:, :, :, :), Val(:), 2) == (1ˢᵗ, 1ˢᵗ, 2ⁿᵈ, 2ⁿᵈ)
end

@testset "Coupling" begin
    tc = Coupling(:, DID, (2,))
    @test tc == Coupling(1, :, DID, (2,)) == Coupling(Index{DID}(:, 2))
    @test id(tc) == tc.pattern
    @test length(tc) == length(typeof(tc)) == 1
    @test eltype(tc) == eltype(typeof(tc)) == typeof(tc)
    @test collect(tc) == [tc]
    @test rank(tc) == rank(typeof(tc)) == 1
    @test tc * tc == tc ⊗ tc
    @test string(tc) == "Index{DID}(:, 2)"
    @test latexstring(tc) == "d^{}_{:,\\,2}"

    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>DFock(2))
    tc₁ = Coupling(1.5, :, DID, (1, 2))
    tc₂ = Coupling(2.0, @pattern(Index{DID}(:, a), Index{DID}(:, b); constraint=a<b))
    ex = expand(tc₁, Val(:), bond, hilbert)
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CoordinatedIndex{Index{DID{Int}, Int}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(1.5, ID(
            CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CoordinatedIndex(Index{DID}(1, 2), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    ]
    ex = expand(tc₂, Val(:), bond, hilbert)
    @test eltype(ex) == eltype(typeof(ex)) == Operator{Float64, NTuple{2, CoordinatedIndex{Index{DID{Int}, Int}, SVector{2, Float64}}}}
    @test collect(ex) == [
        Operator(2.0, ID(
            CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CoordinatedIndex(Index{DID}(1, 2), SVector(0.0, 0.0), SVector(0.0, 0.0))
            )
        )
    ]

    tc = tc₁*tc₂
    @test string(tc) == "3.0 [Index{DID}(:, 1) Index{DID}(:, 2)] ⊗ ∑[Index{DID}(:, a) Index{DID}(:, b)](a < b)"
    @test latexstring(tc) == "3.0\\,d^{}_{:,\\,1} d^{}_{:,\\,2} \\cdot \\sum_{a < b} d^{}_{:,\\,a} d^{}_{:,\\,b}"
end

@testset "MatrixCoupling" begin
    component = Component([1, 2], [2, 1], [-1 0; 0 1])
    @test length(component) == 2
    @test component[1] == (1, 2, -1)
    @test component[2] == (2, 1, +1)

    mc = MatrixCoupling(:, DID, component)
    @test eltype(typeof(mc)) == Coupling{Int64, Pattern{Tuple{Colon, Colon}, InternalPattern{(2,), Tuple{DID{Int}, DID{Int}}, 1, Tuple{AllEqual{()}}}}}
    @test mc[1] == Coupling(-1, Index{DID}(:, 1), Index{DID}(:, 2))
    @test mc[2] == Coupling(+1, Index{DID}(:, 2), Index{DID}(:, 1))
    @test mc^2 == mc*mc
    @test mc/2 == mc*(1/2)
    @test mc//2 == mc*(1//2)
    @test -mc  == (-1)*mc

    another = MatrixCoupling((1ˢᵗ, 2ⁿᵈ), DID, Component([:], [:], hcat(2.0)))
    @test another[1] == Coupling(2.0, Index{DID}(1ˢᵗ, :), Index{DID}(2ⁿᵈ, :))

    mcp = 2 * mc * another
    @test mcp == MatrixCouplingProd(mc, another) * 2
    @test eltype(mcp) == Coupling{Float64, Pattern{Tuple{Colon, Colon, Ordinal, Ordinal}, InternalPattern{(2, 2), Tuple{DID{Int}, DID{Int}, DID{Colon}, DID{Colon}}, 2, Tuple{AllEqual{()}, AllEqual{(:nambu,)}}}}}
    @test mcp[1] == 2 * mc[1] * another[1]
    @test mcp[2] == 2 * mc[2] * another[1]
    @test mc*2 == 2*mc == MatrixCouplingProd(2, mc)
    @test mcp*2 == 2*mcp == MatrixCouplingProd(4, mc, another)
    @test mcp*mc == MatrixCouplingProd(2, mc, another, mc)
    @test mc*mcp == MatrixCouplingProd(2, mc, mc, another)
    @test mcp*mcp == MatrixCouplingProd(4, mc, another, mc, another)
    @test mcp^2 == mcp*mcp
    @test mcp/4 == mcp*(1/4) == MatrixCouplingProd(1/2, mc, another)
    @test mcp//4 == mcp*(1//4) == MatrixCouplingProd(1//2, mc, another)
    @test -mcp == (-1)*mcp

    mc₁ = MatrixCoupling((1ˢᵗ, 2ⁿᵈ), DID, Component([1, 2], [2, 1], [0 1; 1 0]))
    mc₂ = MatrixCoupling((2ⁿᵈ, 1ˢᵗ), DID, Component([1, 2], [2, 1], [0 1im; -1im 0]))
    mcs = mc₁ + mc₂
    @test mcs == MatrixCouplingSum(mc₁, mc₂)
    @test eltype(mcs) == Coupling{Complex{Int64}, Pattern{Tuple{Ordinal, Ordinal}, InternalPattern{(2,), Tuple{DID{Int}, DID{Int}}, 1, Tuple{AllEqual{()}}}}}
    @test collect(mcs) == [
        Coupling(Index{DID}(1ˢᵗ, 2), Index{DID}(2ⁿᵈ, 2)),
        Coupling(Index{DID}(1ˢᵗ, 1), Index{DID}(2ⁿᵈ, 1)),
        Coupling(-1im, Index{DID}(2ⁿᵈ, 2), Index{DID}(1ˢᵗ, 2)),
        Coupling(1im, Index{DID}(2ⁿᵈ, 1), Index{DID}(1ˢᵗ, 1))
    ]
    @test mcs*2 == 2*mcs == MatrixCouplingSum(2*mc₁, 2*mc₂)
    @test mcs*mc₁ == MatrixCouplingSum(mc₁*mc₁, mc₂*mc₁)
    @test mc₂*mcs == MatrixCouplingSum(mc₂*mc₁, mc₂*mc₂)
    @test mcp*mcs == MatrixCouplingSum(2*mc*another*mc₁, 2*mc*another*mc₂)
    @test mcs*mcp == MatrixCouplingSum(2*mc₁*mc*another, 2*mc₂*mc*another)
    @test mcs^2 == mcs*mcs
    @test mcs/2 == mcs*(1/2) == MatrixCouplingSum(mc₁/2, mc₂/2)
    @test mcs//2 == mcs*(1//2) == MatrixCouplingSum(mc₁//2, mc₂//2)
    @test -mcs == (-1)*mcs == MatrixCouplingSum(-mc₁, -mc₂)
    @test mcs+mc₁ == MatrixCouplingSum(mc₁, mc₂, mc₁)
    @test mc₁+mcs == MatrixCouplingSum(mc₁, mc₁, mc₂)
    @test mcs+mcs == MatrixCouplingSum(mc₁, mc₂, mc₁, mc₂)

    mcs₂ = mc₁ - mc₂
    @test mcs₂ == MatrixCouplingSum(mc₁, -mc₂)
    @test mcs - mc₁ == MatrixCouplingSum(mc₁, mc₂, -mc₁)     
    @test mc₁ - mcs₂ == MatrixCouplingSum(mc₁, -mc₁, mc₂)
    @test mcs - mcs₂ == MatrixCouplingSum(mc₁, mc₂, -mc₁, mc₂)
end

@testset "TermFunction" begin
    bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))

    ta = TermAmplitude(nothing)
    @test ta(bond) == 1
    @test valtype(ta, bond) == valtype(typeof(ta), typeof(bond)) == Int

    ta = TermAmplitude(bond::Bond->bond.kind+3.0)
    @test ta(bond) == 4.0
    @test valtype(ta, bond) == valtype(typeof(ta), typeof(bond)) == Float64

    tcs = Coupling(1.0, (1ˢᵗ, 2ⁿᵈ), DID, (1, 1)) + Coupling(2.0, (1ˢᵗ, 2ⁿᵈ), DID, (2, 2))
    termcouplings = TermCoupling(tcs)
    @test termcouplings == deepcopy(TermCoupling(tcs))
    @test isequal(termcouplings, deepcopy(TermCoupling(tcs)))
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == eltype(typeof(tcs))
    @test termcouplings(bond) == tcs

    bond₁ = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
    bond₂ = Bond(2, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))

    fx = bond::Bond -> bond.kind==1 ? Coupling(1.0, (1ˢᵗ, 2ⁿᵈ), DID, (1, 1)) : Coupling(1.0, (1ˢᵗ, 2ⁿᵈ), DID, (2, 2))
    termcouplings = TermCoupling(fx)
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typejoin(typeof(fx(bond₁)), typeof(fx(bond₂)))
    @test termcouplings(bond₁) == fx(bond₁)
    @test termcouplings(bond₂) == fx(bond₂)
end

@testset "Term" begin
    term = Term{:Mu}(:μ, 1.5, 0, bond->iseven(bond[1].site) ? Coupling(1.0, (1ˢᵗ, 1ˢᵗ), DID, (2, 2)) : Coupling(1.0, (1ˢᵗ, 1ˢᵗ), DID, (1, 1)), true; amplitude=bond->3)
    @test term == deepcopy(term)
    @test isequal(term, deepcopy(term))
    @test term|>kind == term|>typeof|>kind == :Mu
    @test term|>id == term|>typeof|>id == :μ
    @test term|>valtype == term|>typeof|>valtype == Float
    @test term|>rank == term|>typeof|>rank == 2

    p₁, p₂ = Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (1.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(DFock(2), 2)
    @test string(term, Bond(p₁), hilbert) == "4.5 Index{DID}(1ˢᵗ, 1) Index{DID}(1ˢᵗ, 1)"
    @test string(term, Bond(p₂), hilbert) == "4.5 Index{DID}(1ˢᵗ, 2) Index{DID}(1ˢᵗ, 2)"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test update!(term, μ=4.25) == replace(term, value=4.25)
    @test term.value == 4.25

    term = Term{:Mu}(:μ, 1.5, 0, Coupling(1.0, :, DID, (2, 1)), true; amplitude=bond->3, ismodulatable=false)
    bond = Bond(Point(1, (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(DFock(2))
    @test string(term, bond, hilbert) == "4.5 Index{DID}(:, 2) Index{DID}(:, 1)"
    operators = Operators(
        Operator(+2.25,
            CoordinatedIndex(Index{DID}(1, 2), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.0), SVector(0.0, 0.0))
        )
    )
    @test expand(term, bond, hilbert, half=true) == expand(term, [bond], hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == expand(term, [bond], hilbert, half=false) == operators*2

    term = Term{:Hp}(:t, 1.5, 1, Coupling(1.0, (1ˢᵗ, 2ⁿᵈ), DID, (2, 1)), false; amplitude=bond->3.0)
    bond = Bond(1, Point(2, (1.5, 1.5), (1.0, 1.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(DFock(2), 2)
    @test string(term, bond, hilbert) == "4.5 Index{DID}(1ˢᵗ, 2) Index{DID}(2ⁿᵈ, 1) + h.c."
    operators = Operators(
        Operator(4.5,
            CoordinatedIndex(Index{DID}(2, 2), SVector(1.5, 1.5), SVector(1.0, 1.0)),
            CoordinatedIndex(Index{DID}(1, 1), SVector(0.5, 0.5), SVector(0.0, 0.0))
        )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Metric" begin
    m = OperatorUnitToTuple((:site, :nambu))
    @test m == OperatorUnitToTuple(:site, :nambu)
    @test isequal(m, OperatorUnitToTuple(:site, :nambu))
    @test keys(m) == keys(typeof(m)) == (:site, :nambu)
    @test OperatorUnitToTuple(Index{DID{Int}, Int}) == OperatorUnitToTuple(:site, :nambu)
    @test OperatorUnitToTuple(CompositeIndex{Index{DID{Int}, Int}}) == OperatorUnitToTuple(:site, :nambu)
    @test OperatorUnitToTuple(Hilbert{DFock}) == OperatorUnitToTuple(:site, :nambu)
    @test valtype(typeof(m), Index{DID{Int}, Int}) == Tuple{Int, Int}
    @test valtype(typeof(m), CompositeIndex{Index{DID{Int}, Int}}) == Tuple{Int, Int}

    index = CoordinatedIndex(Index{DID}(4, 1), SVector(0.5, 0.0), SVector(1.0, 0.0))
    @test m(index.index) == (4, 1) == m(index)
end

@testset "Table" begin
    @test contentnames(Table) == (:by, :contents)

    by = OperatorUnitToTuple(:site)

    table = Table([Index{DID}(1, 1), Index{DID}(1, 2)], by)
    @test empty(table) == Table{Index{DID{Int}, Int}}(by)
    @test table[Index{DID}(1, 1)]==1 && table[Index{DID}(1, 2)]==1

    hilbert = Hilbert(site=>DFock(2) for site in [1, 2])
    inds₁ = (Index(1, internal) for internal in DFock(2))|>collect
    inds₂ = (Index(2, internal) for internal in DFock(2))|>collect
    @test Table(hilbert, by) == Table([inds₁; inds₂], by)
    @test Table(hilbert, by) == union(Table(inds₁, by), Table(inds₂, by))

    opt = Operator(1.0im,
        CoordinatedIndex(Index{DID}(1, 2), SVector(0.0, 0.0), SVector(1.0, 0.0)),
        CoordinatedIndex(Index{DID}(1, 1), SVector(0.0, 0.0), SVector(0.0, 0.0))
        )
    @test sequence(opt, table) == (1, 1)
    @test haskey(table, opt.id) == (true, true)

    table = Table(hilbert)
    @test table == Table([inds₁; inds₂])
    @test reset!(empty(table), [inds₁; inds₂]) == table
    @test reset!(empty(table), hilbert) == table

    hilbert = Hilbert(DFock(2), 5)
    table = Table(hilbert)
    @test findall(index->index.site∈2:4 && index.internal.nambu==1, hilbert, table) == [3, 5, 7]
end

@testset "Boundary" begin
    op = Operator(4.5,
        CoordinatedIndex(Index{DID}(1, 2), SVector(0.5, 0.5), SVector(0.0, 0.0)),
        CoordinatedIndex(Index{DID}(2, 1), SVector(1.5, 1.5), SVector(1.0, 1.0))
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
