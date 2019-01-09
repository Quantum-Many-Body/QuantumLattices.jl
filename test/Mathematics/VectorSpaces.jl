using Hamiltonian.Mathematics.VectorSpaces

@testset "DirectVectorSpace" begin
    id1,id2,id3=(1,1),(1,2),(1,3)

    vs=DirectVectorSpace{'T'}(id1,id2)
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
    @test vs==id1⊕id2
    @test vs⊕id3==id1⊕id2⊕id3
    @test id3⊕vs==id3⊕id1⊕id2
    @test (id2⊕id3)⊕vs==id2⊕id3⊕id1⊕id2

    vs=DirectVectorSpace{'F'}(id1,id2)
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
    @test vs==id1⊕id2
    @test vs⊕id3==id1⊕id2⊕id3
    @test id3⊕vs==id3⊕id1⊕id2
    @test (id2⊕id3)⊕vs==id2⊕id3⊕id1⊕id2
end

@testset "DirectOrderedIndices" begin
    foi=DirectOrderedIndices{'F'}(2,2,2)
    @test dimension(foi)==8
    @test foi|>collect==[(1,1,1),(2,1,1),(1,2,1),(2,2,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]
    @test (1,1,1) ∈ foi && (1,2,3) ∉ foi
    for i=1:length(foi)
        @test searchsortedfirst(foi,foi[i])==i
    end
    coi=DirectOrderedIndices{'C'}(2,2,2)
    @test dimension(coi)==8
    @test (1,1,1) ∈ coi && (1,2,3) ∉ coi
    @test coi|>collect==[(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    for i=1:length(coi)
        @test searchsortedfirst(coi,coi[i])==i
    end
end

@testset "TabledOrderedIndices" begin
    dims=(2,2)
    table=[(1,1),(1,2),(2,1),(2,2)]
    toi=TabledOrderedIndices{'T'}(dims,table)
    @test toi==TabledOrderedIndices{2}(DulPermutations,2)
    @test dimension(toi)==4
    @test toi|>collect==table
    @test (1,1) ∈ toi && (1,3) ∉ toi
    for i=1:length(toi)
        @test searchsortedfirst(toi,toi[i])==i
    end

    table=[(1,2),(2,1),(2,2),(1,1)]
    toi=TabledOrderedIndices{'F'}(dims,table)
    @test dimension(toi)==4
    @test toi|>collect==table
    @test (1,1) ∈ toi && (1,3) ∉ toi
    for i=1:length(toi)
        @test searchsortedfirst(toi,toi[i])==i
    end
end
