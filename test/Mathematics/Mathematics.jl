@time @safetestset "Mathematics" begin
    @time @safetestset "Combinatorics" begin include("Combinatorics.jl") end
    @time @safetestset "VectorSpaces" begin include("VectorSpaces.jl") end
    @time @safetestset "AlgebraOverFields" begin include("AlgebraOverFields.jl") end
    @time @safetestset "QuantumNumbers" begin include("QuantumNumbers.jl") end
end
