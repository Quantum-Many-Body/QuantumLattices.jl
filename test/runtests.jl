using Test
using SafeTestsets

@testset "all" begin
    include("Prerequisites/Prerequisites.jl")
    include("Mathematics/Mathematics.jl")
    include("Essentials/Essentials.jl")
end
