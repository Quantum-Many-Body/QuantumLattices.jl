@time @safetestset "Essentials" begin
    @time @safetestset "QuantumOperators" begin include("QuantumOperators.jl") end
    @time @safetestset "QuantumNumbers" begin include("QuantumNumbers.jl") end
    @time @safetestset "Spatials" begin include("Spatials.jl") end
    @time @safetestset "DegreesOfFreedom" begin include("DegreesOfFreedom.jl") end
    @time @safetestset "QuantumSystems" begin include("QuantumSystems.jl") end
    @time @safetestset "Frameworks" begin include("Frameworks.jl") end
end
