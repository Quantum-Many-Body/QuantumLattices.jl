using Test
using QuantumLattices.Prerequisites.Trees

@tree struct ATree end

@testset "AbstractTree" begin
    tree=ATree{String,Int}()

    @test tree|>eltype==Pair{String,Int}
    @test tree|>keytype==String
    @test tree|>valtype==Int
    @test tree|>root==nothing

    push!(tree,"L0",0)
    push!(tree,"L0","L1-1",2)
    push!(tree,"L0","L1-2",3)
    push!(tree,"L1-1","L2-1",4)
    push!(tree,"L1-1","L2-2",5)
    push!(tree,"L1-2","L2-3",6)
    push!(tree,"L1-2","L2-4",7)
    tree["L0"]=1

    @test tree|>root=="L0"
    @test tree|>length==7
    @test keys(tree,treedepth)|>collect==["L0","L1-1","L2-1","L2-2","L1-2","L2-3","L2-4"]
    @test keys(tree,treewidth)|>collect==["L0","L1-1","L1-2","L2-1","L2-2","L2-3","L2-4"]
    @test values(tree,treedepth)|>collect==[1,2,4,5,3,6,7]
    @test values(tree,treewidth)|>collect==[1,2,3,4,5,6,7]
    @test pairs(tree,treedepth)|>collect==[("L0",1),("L1-1",2),("L2-1",4),("L2-2",5),("L1-2",3),("L2-3",6),("L2-4",7)]
    @test pairs(tree,treewidth)|>collect==[("L0",1),("L1-1",2),("L1-2",3),("L2-1",4),("L2-2",5),("L2-3",6),("L2-4",7)]
    @test [isleaf(tree,node) for node in keys(tree,treewidth)]==[false,false,false,true,true,true,true]
    @test [level(tree,node) for node in keys(tree,treewidth)]==[1,2,2,3,3,3,3]
    @test [ancestor(tree,"L2-1",i) for i=0:2]==["L2-1","L1-1","L0"]
    @test [ancestor(tree,"L2-3",i) for i=0:2]==["L2-3","L1-2","L0"]
    @test descendants(tree,"L0",0)==["L0"]
    @test descendants(tree,"L0",1)==["L1-1","L1-2"]
    @test descendants(tree,"L0",2)==["L2-1","L2-2","L2-3","L2-4"]
    @test siblings(tree,"L0")==[]
    @test siblings(tree,"L1-1")==["L1-2"]
    @test siblings(tree,"L2-1")==["L2-2"]
    @test leaves(tree)==["L2-1","L2-2","L2-3","L2-4"]

    backup=deepcopy(tree)

    sub=subtree(tree,"L0")
    empty!(tree)
    append!(tree,sub)
    @test tree==backup
    @test isequal(tree,backup)

    sub=subtree(tree,"L1-1")
    delete!(tree,"L1-1")
    append!(tree,"L0",sub)
    move!(tree,"L1-2","L0")
    @test tree==backup
    @test isequal(tree,backup)
end

@testset "SimpleTree" begin
    tree=SimpleTree{String,Int}()
    @test tree|>eltype==Pair{String,Int}
    @test tree|>keytype==String
    @test tree|>valtype==Int
    @test tree|>root==nothing
end
