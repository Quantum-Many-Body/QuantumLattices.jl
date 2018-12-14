@testset "essentials" begin
    @testset "Spatial" begin include("Spatial.jl") end
    @testset "DegreeOfFreedom" begin include("DegreeOfFreedom.jl") end
end
