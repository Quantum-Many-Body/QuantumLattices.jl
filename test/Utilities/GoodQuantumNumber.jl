using Test
using Hamiltonian.Utilities.GoodQuantumNumber

@testset "SQN" begin
    @test SQN|>fieldnames==(:Sz,)
    @test SQN|>periods==(Inf,)
    @test SQN(0.0)==SQN(0.0)
    @test SQN(1.0)|>zero==SQN(0.0)
    @test string(SQN(0.0))=="SQN(0.0)"
    @test length(SQN(0.0))==1
    @test SQN(0.0).Sz==0.0
    @test regularization(SQN,[1.5])==[1.5]
    @test SQN(0.5)+SQN(-0.5)==SQN(0.0)
    @test SQN(0.5)-SQN(-0.5)==SQN(1.0)
    @test 2*SQN(0.5)==SQN(0.5)*2==SQN(1.0)
    @test replace(SQN(0.5),Sz=1.5)==SQN(1.5)
end

@testset "PQN" begin
    @test PQN|>fieldnames==(:N,)
    @test PQN|>periods==(Inf,)
    @test PQN(0.0)==PQN(0.0)
    @test PQN(1.0)|>zero==PQN(0.0)
    @test string(PQN(0.0))=="PQN(0.0)"
    @test length(PQN(0.0))==1
    @test PQN(0.0).N==0.0
    @test regularization(PQN,[2.0])==[2.0]
    @test PQN(1.0)+PQN(2.0)==PQN(3.0)
    @test PQN(1.0)-PQN(2.0)==PQN(-1.0)
    @test PQN(1.0)*2==2*PQN(1.0)==PQN(2.0)
    @test replace(PQN(0.0),N=4.0)==PQN(4.0)
end

@testset "SPQN" begin
    @test SPQN|>fieldnames==(:N,:Sz)
    @test SPQN|>periods==(Inf,Inf)
    @test SPQN(0.0,0.0)==SPQN(0.0,0.0)
    @test SPQN(1.0,0.5)|>zero==SPQN(0.0,0.0)
    @test string(SPQN(0.0,0.0))=="SPQN(0.0,0.0)"
    @test length(SPQN(0.0,0.0))==2
    @test SPQN(0.0,0.0).N==0.0
    @test SPQN(0.0,0.0).Sz==0.0
    @test regularization(SPQN,[2.0,1.5])==[2.0,1.5]
    @test SPQN(1.0,0.5)+SPQN(1.0,-0.5)==SPQN(2.0,0.0)
    @test SPQN(1.0,0.5)-SPQN(1.0,-0.5)==SPQN(0.0,1.0)
    @test SPQN(1.0,0.5)*2==2*SPQN(1.0,0.5)==SPQN(2.0,1.0)
    @test replace(SPQN(1.0,0.5),N=3.0)==SPQN(3.0,0.5)
    @test replace(SPQN(1.0,0.5),Sz=-0.5)==SPQN(1.0,-0.5)
    @test replace(SPQN(1.0,0.5),N=2.0,Sz=1.0)==SPQN(2.0,1.0)
end

@testset "Z2QN" begin
    @test Z2QN|>fieldnames==(:N,)
    @test Z2QN|>periods==(2,)
    @test Z2QN(0.0)==Z2QN(0.0)
    @test Z2QN(1.0)|>zero==Z2QN(0.0)
    @test string(Z2QN(0.0))=="Z2QN(0.0)"
    @test length(Z2QN(0.0))==1
    @test Z2QN(0.0).N==0.0
    @test regularization(Z2QN,[2.0])==[0.0]
    @test Z2QN(1.0)+Z2QN(1.0)==Z2QN(0.0)
    @test Z2QN(1.0)-Z2QN(1.0)==Z2QN(0.0)
    @test Z2QN(1.0)*2==2*Z2QN(1.0)==Z2QN(0.0)
    @test replace(Z2QN(1.0),N=0.0)==Z2QN(0.0)
end

@testset "quantumnumber" begin
    qndict=Dict(SPQN(1.0,-0.5)=>1,SPQN(1.0,0.5)=>2)
    @test qndict[SPQN(1.0,-0.5)]==1
    @test qndict[SPQN(1.0,0.5)]==2
    @test ⊗(SPQN,PQN(1.0),SQN(-0.5))==SPQN(1.0,-0.5)
end

@testset "quantumnumbers" begin
    println(QuantumNumbers('C',SQN,[[1.0] [2.0]],[1,1],qnscounts))
    println(QuantumNumbers('C',[SQN(1.0),SQN(2.0)],[1,1],qnscounts))
    println(QuantumNumbers(SQN(1.0),4))
end
