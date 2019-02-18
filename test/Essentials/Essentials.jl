@testset "Essentials" begin
    @testset "Spatials" begin include("Spatials.jl") end
    @testset "DegreesOfFreedom" begin include("DegreesOfFreedom.jl") end
    @testset "Terms" begin include("Terms.jl") end
    @testset "FockPackage" begin include("FockPackage.jl") end
    @testset "SpinPackage" begin include("SpinPackage.jl") end
end
