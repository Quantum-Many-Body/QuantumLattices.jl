@time @safetestset "Essentials" begin
    @time @safetestset "QuantumAlgebras" begin include("QuantumAlgebras.jl") end
    @time @safetestset "QuantumNumbers" begin include("QuantumNumbers.jl") end
    @time @safetestset "Spatials" begin include("Spatials.jl") end
    @time @safetestset "DegreesOfFreedom" begin include("DegreesOfFreedom.jl") end
    @time @safetestset "Terms" begin include("Terms.jl") end
    @time @safetestset "Frameworks" begin include("Frameworks.jl") end
    @time @safetestset "QuantumSystems" begin include("QuantumSystems.jl") end
end
