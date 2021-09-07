@time @safetestset "Essentials" begin
    @time @safetestset "Spatials" begin include("Spatials.jl") end
    @time @safetestset "DegreesOfFreedom" begin include("DegreesOfFreedom.jl") end
    @time @safetestset "Terms" begin include("Terms.jl") end
    @time @safetestset "Frameworks" begin include("Frameworks.jl") end
    @time @safetestset "FockPackage" begin include("FockPackage.jl") end
    @time @safetestset "SpinPackage" begin include("SpinPackage.jl") end
    @time @safetestset "PhononPackage" begin include("PhononPackage.jl") end
end
