using Test
using QuantumLattices.Mathematics.VectorSpaces
using QuantumLattices.Mathematics.Combinatorics: Combinations
using QuantumLattices.Interfaces: dimension,⊕,rank,dims,inds

@testset "SimpleVectorSpace" begin
    id1,id2,id3=(1,1),(1,2),(1,3)

    vs=SimpleVectorSpace{'T'}(id1,id2)
    @test vs==deepcopy(vs)
    @test isequal(vs,deepcopy(vs))
    @test vs|>size==(2,)
    @test vs|>dimension==2
    @test vs|>collect==[id1,id2]
    @test vs[1]==id1 && vs[2]==id2
    @test searchsortedfirst(vs,id1)==1 && searchsortedfirst(vs,id2)==2
    @test findfirst(id1,vs)==1 && findfirst(id2,vs)==2
    @test findfirst((id1,id2),vs)==(1,2)
    @test id1 ∈ vs && id2 ∈ vs && id3 ∉ vs
    @test HasTable(typeof(vs))==HasTable(true)
    @test TableSorted(typeof(vs))==TableSorted(true)
    @test IsMultiIndexable(typeof(vs))==IsMultiIndexable(false)
    @test vs==id1⊕id2
    @test vs⊕id3==id1⊕id2⊕id3
    @test id3⊕vs==id3⊕id1⊕id2
    @test (id2⊕id3)⊕vs==id2⊕id3⊕id1⊕id2

    vs=SimpleVectorSpace{'F'}(id1,id2)
    @test vs==deepcopy(vs)
    @test isequal(vs,deepcopy(vs))
    @test vs|>size==(2,)
    @test vs|>dimension==2
    @test vs|>collect==[id1,id2]
    @test vs[1]==id1 && vs[2]==id2
    @test searchsortedfirst(vs,id1)==1 && searchsortedfirst(vs,id2)==2
    @test findfirst(id1,vs)==1 && findfirst(id2,vs)==2
    @test findfirst((id1,id2),vs)==(1,2)
    @test id1 ∈ vs && id2 ∈ vs && id3 ∉ vs
    @test HasTable(typeof(vs))==HasTable(true)
    @test TableSorted(typeof(vs))==TableSorted(false)
    @test IsMultiIndexable(typeof(vs))==IsMultiIndexable(false)
    @test vs==id1⊕id2
    @test vs⊕id3==id1⊕id2⊕id3
    @test id3⊕vs==id3⊕id1⊕id2
    @test (id2⊕id3)⊕vs==id2⊕id3⊕id1⊕id2
end

@testset "SimpleIndices" begin
    foi=SimpleIndices{'F'}(2,2,2)
    @test HasTable(typeof(foi))==HasTable(false)
    @test IsMultiIndexable(typeof(foi))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(foi))==MultiIndexOrderStyle('F')
    @test dimension(foi)==8
    @test dims(foi)==(2,2,2)
    @test rank(typeof(foi))==3
    @test inds((1,1,1),foi)==(1,1,1)
    @test Tuple((1,1,1),foi)==(1,1,1)
    @test foi|>collect==[(1,1,1),(2,1,1),(1,2,1),(2,2,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]
    @test (1,1,1) ∈ foi && (1,2,3) ∉ foi
    for (i,finds) in enumerate(foi)
        @test findfirst(finds,foi)==i
        @test searchsortedfirst(foi,foi[i])==i
    end
    coi=SimpleIndices{'C'}(2,2,2)
    @test HasTable(typeof(coi))==HasTable(false)
    @test IsMultiIndexable(typeof(coi))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(coi))==MultiIndexOrderStyle('C')
    @test dimension(coi)==8
    @test dims(coi)==(2,2,2)
    @test rank(typeof(coi))==3
    @test inds((1,1,1),coi)==(1,1,1)
    @test Tuple((1,1,1),coi)==(1,1,1)
    @test (1,1,1) ∈ coi && (1,2,3) ∉ coi
    @test coi|>collect==[(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    for (i,cinds) in enumerate(coi)
        @test findfirst(cinds,coi)==i
        @test searchsortedfirst(coi,coi[i])==i
    end
end

@testset "TabledIndices" begin
    dims=(2,2)
    table=[(1,1),(1,2),(2,1),(2,2)]
    toi=TabledIndices{'T'}(dims,table)
    @test HasTable(typeof(toi))==HasTable(true)
    @test TableSorted(typeof(toi))==TableSorted(true)
    @test toi==TabledIndices{2}(DulPermutations,2)
    @test dimension(toi)==4
    @test toi|>collect==table
    @test (1,1) ∈ toi && (1,3) ∉ toi
    for i=1:length(toi)
        @test searchsortedfirst(toi,toi[i])==i
    end

    table=[(1,2),(2,1),(2,2),(1,1)]
    toi=TabledIndices{'F'}(dims,table)
    @test HasTable(typeof(toi))==HasTable(true)
    @test TableSorted(typeof(toi))==TableSorted(false)
    @test dimension(toi)==4
    @test toi|>collect==table
    @test (1,1) ∈ toi && (1,3) ∉ toi
    for i=1:length(toi)
        @test searchsortedfirst(toi,toi[i])==i
    end
end

struct VSZNamedVectorSpace{NS,BS<:Tuple,VS<:Tuple{Vararg{Vector}}} <: NamedVectorSpace{:zip,NS,BS,VS}
    contents::VS
end
@generated function VSZNamedVectorSpace{NS}(contents::Vector...) where NS
    @assert length(NS)==length(contents) && isa(NS,Tuple{Vararg{Symbol}})
    BS=Expr(:curly,:Tuple,[contents[i]|>eltype for i=1:length(NS)]...)
    return quote
        @assert mapreduce(length,==,contents)
        VSZNamedVectorSpace{NS,$BS,typeof(contents)}(contents)
    end
end

struct VSPNamedVectorSpace{NS,BS<:Tuple,VS<:Tuple{Vararg{Vector}}} <: NamedVectorSpace{:product,NS,BS,VS}
    contents::VS
end
@generated function VSPNamedVectorSpace{NS}(contents::Vector...) where NS
    @assert length(NS)==length(contents) && isa(NS,Tuple{Vararg{Symbol}})
    BS=Expr(:curly,:Tuple,[contents[i]|>eltype for i=1:length(NS)]...)
    return :(VSPNamedVectorSpace{NS,$BS,typeof(contents)}(contents))
end

@testset "NamedVectorSpace" begin
    @test IsMultiIndexable(NamedVectorSpace)==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(NamedVectorSpace)==MultiIndexOrderStyle('C')

    nvs=VSZNamedVectorSpace{(:t,:U)}([1,2],[8.0,9.0])
    @test nvs|>keys==nvs|>typeof|>keys==(:t,:U)
    @test nvs|>values==([1,2],[8.0,9.0])
    @test nvs|>pairs|>collect==[:t=>[1,2],:U=>[8.0,9.0]]
    @test eltype(nvs,1)==eltype(nvs|>typeof,1)==Int
    @test eltype(nvs,2)==eltype(nvs|>typeof,2)==Float64
    @test nvs|>typeof|>rank==1
    @test dims(nvs)==(2,)
    elements=[(t=1,U=8.0),(t=2,U=9.0)]
    for i=1:dimension(nvs)
        @test NamedTuple(inds(elements[i],nvs),nvs)==elements[i]
    end
    @test nvs|>collect==elements

    nvs=VSPNamedVectorSpace{(:t,:U)}([1.0,2.0],[8.0,9.0])
    @test nvs|>typeof|>rank==2
    @test dims(nvs)==(2,2)
    elements=[(t=1.0,U=8.0),(t=1.0,U=9.0),(t=2.0,U=8.0),(t=2.0,U=9.0)]
    for i=1:dimension(nvs)
        @test NamedTuple(inds(elements[i],nvs),nvs)==elements[i]
    end
    @test nvs|>collect==elements
end
