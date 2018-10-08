using Test
using Hamiltonian.Utilities.GoodQuantumNumber

@testset "quantumnumber" begin
    @test SPQN|>zero==SPQN(0.0,0.0)
    @test SPQN|>fieldnames==(:N,:Sz)
    @test SPQN|>periods==(nothing,nothing)

    qn=SPQN(0.0,0.0)
    @test qn==qn
    @test string(qn)=="SPQN(0.0,0.0)"
    @test length(qn)==2
    @test qn.N==0.0
    @test qn.Sz==0.0

    qndict=Dict(SPQN(1.0,-0.5)=>1,SPQN(1.0,0.5)=>2)
    a1,a2=SPQN(1.0,-0.5),SPQN(1.0,0.5)
    @test qndict[a1]==1
    @test qndict[a2]==2
end
