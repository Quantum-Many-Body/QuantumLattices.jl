@testset "Essentials" begin
    @testset "Spatials" begin include("Spatials.jl") end
    @testset "DegreesOfFreedom" begin include("DegreesOfFreedom.jl") end
    include("FockPackage/FockPackage.jl")
end
