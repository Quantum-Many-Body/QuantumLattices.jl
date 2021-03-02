@time @safetestset "Prerequisites" begin
    @time @safetestset "Tools" begin include("Tools.jl") end
    @time @safetestset "TypeTraits" begin include("TypeTraits.jl") end
    @time @safetestset "CompositeStructures" begin include("CompositeStructures.jl") end
    @time @safetestset "SimpleTrees" begin include("SimpleTrees.jl") end
    @time @safetestset "NamedVectors" begin include("NamedVectors.jl") end
end
