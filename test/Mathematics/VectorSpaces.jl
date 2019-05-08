using Test
using QuantumLattices.Mathematics.VectorSpaces
using QuantumLattices.Mathematics.Combinatorics: Combinations
using QuantumLattices.Interfaces: dimension,⊕,rank,dims,inds,degree

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

@testset "GradedTables" begin
    com=Combinations
    t=GradedTables(com,2,Val((0,1,2)))
    @test keytype(t,1)==keytype(t|>typeof,1)==Int
    @test keytype(t,2)==keytype(t|>typeof,2)==Int
    @test keytype(t,3)==keytype(t|>typeof,3)==Int
    @test valtype(t,1)==valtype(t|>typeof,1)==TabledIndices{'T',0}
    @test valtype(t,2)==valtype(t|>typeof,2)==TabledIndices{'T',1}
    @test valtype(t,3)==valtype(t|>typeof,3)==TabledIndices{'T',2}
    @test rank(t)==rank(typeof(t))==3
    @test t[1]==TabledIndices{0}(com,2)
    @test t[2]==TabledIndices{1}(com,2)
    @test t[3]==TabledIndices{2}(com,2)
end

struct SimpleGradedVectorSpace{T<:GradedTables{Char,SimpleVectorSpace{'F',Int,4}}} <: GradedVectorSpace{Char,Int,SimpleVectorSpace{'F',Int,4},T}
    tables::T
end

@testset "GradedVectorSpace" begin
    v1=SimpleVectorSpace{'F'}(2,3,4,1)
    v2=SimpleVectorSpace{'F'}(8,5,6,7)
    vs=SimpleGradedVectorSpace(GradedTables((v1,v2),('a','b')))
    @test keys(vs)==('a','b')
    @test values(vs)==(v1,v2)
    @test pairs(vs)==vs.tables
    @test keytype(vs,1)==keytype(typeof(vs),1)==Char
    @test keytype(vs,2)==keytype(typeof(vs),2)==Char
    @test valtype(vs,1)==valtype(typeof(vs),1)==SimpleVectorSpace{'F',Int,4}
    @test valtype(vs,2)==valtype(typeof(vs),2)==SimpleVectorSpace{'F',Int,4}
    @test eltype(vs,1)==eltype(typeof(vs),1)==Tuple{Char,Int}
    @test eltype(vs,2)==eltype(typeof(vs),2)==Tuple{Char,Int}
    @test rank(vs)==rank(typeof(vs))==2
    @test degree('a',vs)==1
    @test degree('b',vs)==2
    @test dimension(vs)==dimension(vs,('a','b'))==8
    @test dimension(vs,'a')==dimension(vs,'b')==4
    elements=[('a',2),('a',3),('a',4),('a',1),('b',8),('b',5),('b',6),('b',7)]
    for (i,e) in enumerate(vs)
        @test e==elements[i]
        @test searchsortedfirst(vs,e)==i
        @test findfirst(e,vs)==i
    end
    @test vs|>collect==elements
    for i=1:4
        @test vs[('a',i)]==v1[i]
        @test vs[('b',i)]==v2[i]
    end
end
