using Test
using QuantumLattices.Prerequisites.SimpleTrees
using QuantumLattices.Prerequisites.Traits: contentnames

struct Tree{N, D} <: AbstractSimpleTree{N, D}
    TREECORE::SimpleTreeCore{N, D}
end
Tree{N, D}() where {N, D} = Tree(SimpleTreeCore{N, D}())

@testset "AbstractSimpleTree" begin
    tree = Tree{String, Int}()

    @test contentnames(typeof(tree)) == (:TREECORE,)
    @test eltype(tree) == Pair{String, Int}
    @test keytype(tree) == String
    @test valtype(tree) == Int
    @test root(tree) === nothing

    push!(tree, "L0", 0)
    push!(tree, "L0", "L1-1", 2)
    push!(tree, "L0", "L1-2", 3)
    push!(tree, "L1-1", "L2-1", 4)
    push!(tree, "L1-1", "L2-2", 5)
    push!(tree, "L1-2", "L2-3", 6)
    push!(tree, "L1-2", "L2-4", 7)
    tree["L0"] = 1

    @test root(tree) == "L0"
    @test length(tree) == 7
    @test collect(keys(tree, simpletreedepth)) == ["L0", "L1-1", "L2-1", "L2-2", "L1-2", "L2-3", "L2-4"]
    @test collect(keys(tree, simpletreewidth)) == ["L0", "L1-1", "L1-2", "L2-1", "L2-2", "L2-3", "L2-4"]
    @test collect(values(tree, simpletreedepth)) == [1, 2, 4, 5, 3, 6, 7]
    @test collect(values(tree, simpletreewidth)) == [1, 2, 3, 4, 5, 6, 7]
    @test collect(pairs(tree, simpletreedepth)) == [("L0", 1), ("L1-1", 2), ("L2-1", 4), ("L2-2", 5), ("L1-2", 3), ("L2-3", 6), ("L2-4", 7)]
    @test collect(pairs(tree, simpletreewidth)) == [("L0", 1), ("L1-1", 2), ("L1-2", 3), ("L2-1", 4), ("L2-2", 5), ("L2-3", 6), ("L2-4", 7)]
    @test [isleaf(tree, node) for node in keys(tree, simpletreewidth)] == [false, false, false, true, true, true, true]
    @test [level(tree, node) for node in keys(tree, simpletreewidth)] == [1, 2, 2, 3, 3, 3, 3]
    @test [ancestor(tree, "L2-1", i) for i = 0:2] == ["L2-1", "L1-1", "L0"]
    @test [ancestor(tree, "L2-3", i) for i = 0:2] == ["L2-3", "L1-2", "L0"]
    @test descendants(tree, "L0", 0) == ["L0"]
    @test descendants(tree, "L0", 1) == ["L1-1", "L1-2"]
    @test descendants(tree, "L0", 2) == ["L2-1", "L2-2", "L2-3", "L2-4"]
    @test siblings(tree, "L0") == []
    @test siblings(tree, "L1-1") == ["L1-2"]
    @test siblings(tree, "L2-1") == ["L2-2"]
    @test leaves(tree) == ["L2-1", "L2-2", "L2-3", "L2-4"]

    backup = deepcopy(tree)

    sub = subtree(tree, "L0")
    empty!(tree)
    append!(tree, sub)
    @test tree == backup
    @test isequal(tree, backup)

    sub = subtree(tree, "L1-1")
    delete!(tree, "L1-1")
    append!(tree, "L0", sub)
    move!(tree, "L1-2", "L0")
    @test tree == backup
    @test isequal(tree, backup)
end

@testset "SimpleTree" begin
    tree = SimpleTree{String, Int}()
    @test eltype(tree)  == eltype(typeof(tree)) == Pair{String, Int}
    @test keytype(tree) == keytype(typeof(tree)) == String
    @test valtype(tree) == valtype(typeof(tree)) == Int
    @test root(tree) === nothing
end
