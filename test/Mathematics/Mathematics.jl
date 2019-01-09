@testset "Mathematics" begin
    @testset "Combinatorics" begin include("Combinatorics.jl") end
    @testset "VectorSpaces" begin include("VectorSpaces.jl") end
    @testset "AlgebraOverFields" begin include("AlgebraOverFields.jl") end
    @testset "QuantumNumbers" begin include("QuantumNumbers.jl") end
end
