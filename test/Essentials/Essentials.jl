@testset "Essentials" begin
    @testset "Spatial" begin include("Spatial.jl") end
    @testset "DegreeOfFreedom" begin include("DegreeOfFreedom.jl") end
    include("FockPackage/FockPackage.jl")
end
