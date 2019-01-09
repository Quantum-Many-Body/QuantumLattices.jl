@testset "Prerequisites" begin
    @testset "Tools" begin include("Tools.jl") end
    @testset "TypeTraits" begin include("TypeTraits.jl") end
    @testset "Factories" begin include("Factories.jl") end
    @testset "CompositeStructures" begin include("CompositeStructures.jl") end
    @testset "Trees" begin include("Trees.jl") end
    @testset "NamedVectors" begin include("NamedVectors.jl") end
end
