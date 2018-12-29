@testset "Utilities" begin
    @testset "Tools" begin include("Tools.jl") end
    @testset "TypeTrait" begin include("TypeTrait.jl") end
    @testset "Factory" begin include("Factory.jl") end
    @testset "CompositeStructure" begin include("CompositeStructure.jl") end
    @testset "Tree" begin include("Tree.jl") end
    @testset "NamedVector" begin include("NamedVector.jl") end
    @testset "Combinatorics" begin include("Combinatorics.jl") end
    @testset "AlgebraOverField" begin include("AlgebraOverField.jl") end
    @testset "QuantumNumber" begin include("QuantumNumber.jl") end
end
