@time @safetestset "Prerequisites" begin
    @time @safetestset "Tools" begin include("Tools.jl") end
    @time @safetestset "Combinatorics" begin include("Combinatorics.jl") end
    @time @safetestset "Traits" begin include("Traits.jl") end
    @time @safetestset "CompositeStructures" begin include("CompositeStructures.jl") end
    @time @safetestset "SimpleTrees" begin include("SimpleTrees.jl") end
    @time @safetestset "NamedVectors" begin include("NamedVectors.jl") end
    @time @safetestset "VectorSpaces" begin include("VectorSpaces.jl") end
end
