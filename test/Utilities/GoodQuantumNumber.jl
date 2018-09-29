using Test
using Hamiltonian.Utilities.GoodQuantumNumber

@testset "quantumnumber" begin
    a1,a2=SPQN(1.0,-0.5),SPQN(1.0,0.5)
    println(a1)
    println(a2)
end
