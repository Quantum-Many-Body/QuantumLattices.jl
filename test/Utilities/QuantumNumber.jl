using Test
using Hamiltonian.Utilities.QuantumNumber
using DataStructures: OrderedDict
using Printf: @sprintf
using Base.Iterators: reverse

@testset "SQN" begin
    @test SQN|>fieldnames==(:Sz,)
    @test SQN|>periods==(Inf,)
end

@testset "PQN" begin
    @test PQN|>fieldnames==(:N,)
    @test PQN|>periods==(Inf,)
end

@testset "SPQN" begin
    @test SPQN|>fieldnames==(:N,:Sz)
    @test SPQN|>periods==(Inf,Inf)
end

@testset "Z2QN" begin
    @test Z2QN|>fieldnames==(:N,)
    @test Z2QN|>periods==(2,)
end

@testset "SQNS" begin
    @test SQNS(0.5)==QuantumNumbers('C',[SQN(-0.5),SQN(0.5)],[0,1,2],qnsindptr)
end

@testset "PQNS" begin
    @test PQNS(1.0)==QuantumNumbers('C',[PQN(0.0),PQN(1.0)],[0,1,2],qnsindptr)
end

@testset "SzPQNS" begin
    @test SzPQNS(0.5)==QuantumNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,0.5)],[0,1,2],qnsindptr)
end

@testset "SPQNS" begin
    @test SPQNS(0.5)==QuantumNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,-0.5),SPQN(1.0,0.5),SPQN(2.0,0.0)],[0,1,2,3,4],qnsindptr)
end

@testset "Z2QNS" begin
    @test Z2QNS()==QuantumNumbers('C',[Z2QN(0.0),Z2QN(1.0)],[0,1,2],qnsindptr)
end

@quantumnumber "CN" (:N,) (Inf,)
@quantumnumber "Z4" (:Z,) (4,)
@quantumnumber "CNZ4" (:N,:Z) (Inf,4)

@testset "regularize" begin
    @test regularize(CNZ4,[1.5,5.0])==[1.5,1.0]
    @test regularize(CNZ4,[4.0,-1.0])==[4.0,3.0]
    @test regularize(CNZ4,[1.5 4.0; 5.0 -1.0])==[1.5 4.0; 1.0 3.0]
end

@testset "QuantumNumbers" begin
    qn1,qn2=CNZ4(+1.0,+3.0),CNZ4(-1.0,+1.0)
    qns1=QuantumNumbers('U',[qn1,qn2],[0,2,5],qnsindptr)
    qns2=QuantumNumbers('U',[qn1,qn2],[2,3],qnscounts)
    qns3=QuantumNumbers(OrderedDict(qn1=>2,qn2=>3))
    qns4=QuantumNumbers(OrderedDict(qn1=>1:2,qn2=>3:5))
    @test qns1==qns2==qns3==qns4

    qns=QuantumNumbers('U',[qn1,qn2],[0,3,5],qnsindptr)
    @test qns|>string=="QNS(2,5)"
    @test @sprintf("%s",qns)=="QNS(CNZ4(1.0,3.0)=>1:3,CNZ4(-1.0,1.0)=>4:5)"
    @test qns[1]==qn1
    @test qns[2]==qn2
    @test qns[1:2]==QuantumNumbers('G',[qn1,qn2],[0,3,5],qnsindptr)
    @test qns[[2,1]]==QuantumNumbers('G',[qn2,qn1],[0,2,5],qnsindptr)
    @test qns|>collect==[qn1,qn2]
    @test qns|>reverse|>collect==[qn2,qn1]
    @test qns|>keys|>collect==[qn1,qn2]
    @test values(qns,qnsindptr)|>collect==[1:3,4:5]
    @test values(qns,qnscounts)|>collect==[3,2]
    @test pairs(qns,qnsindptr)|>collect==[qn1=>1:3,qn2=>4:5]
    @test pairs(qns,qnscounts)|>collect==[qn1=>3,qn2=>2]
end

@testset "arithmetic" begin
    qn=CNZ4(1.0,3.0)
    @test +qn==CNZ4(+1.0,+3.0)
    @test -qn==CNZ4(-1.0,+1.0)
    @test qn*4==4*qn==CNZ4(4.0,0.0)
    @test qn^3==CNZ4(3.0,1.0)
    @test ⊗(CNZ4,CN(1.0),Z4(3.0))==qn

    qn1,qn2=CNZ4(1.0,2.0),CNZ4(2.0,3.0)
    @test qn1+qn2==CNZ4(3.0,1.0)
    @test qn1-qn2==CNZ4(-1.0,3.0)
    @test ⊗((qn1,qn2),(+1,+1))==CNZ4(+3.0,1.0)
    @test ⊗((qn1,qn2),(+1,-1))==CNZ4(-1.0,3.0)
    @test ⊗((qn1,qn2),(-1,+1))==CNZ4(+1.0,1.0)
    @test ⊗((qn1,qn2),(-1,-1))==CNZ4(-3.0,3.0)
    @test ⊕((qn1,qn2),(+1,+1))==QuantumNumbers('G',[+qn1,+qn2],[0,1,2],qnsindptr)
    @test ⊕((qn1,qn2),(+1,-1))==QuantumNumbers('G',[+qn1,-qn2],[0,1,2],qnsindptr)
    @test ⊕((qn1,qn2),(-1,+1))==QuantumNumbers('G',[-qn1,+qn2],[0,1,2],qnsindptr)
    @test ⊕((qn1,qn2),(-1,-1))==QuantumNumbers('G',[-qn1,-qn2],[0,1,2],qnsindptr)

    qns=QuantumNumbers('U',[qn1,qn2],[0,2,4],qnsindptr)
    @test +qns==QuantumNumbers('U',[+qn1,+qn2],[0,2,4],qnsindptr)
    @test -qns==QuantumNumbers('U',[-qn1,-qn2],[0,2,4],qnsindptr)
    @test qns*3==3*qns==QuantumNumbers('U',[qn1*3,qn2*3],[0,2,4],qnsindptr)
    @test qns^2==QuantumNumbers('G',[qn1+qn1,qn1+qn2,qn1+qn1,qn1+qn2,qn2+qn1,qn2+qn2,qn2+qn1,qn2+qn2],[2,2,2,2,2,2,2,2],qnscounts)

    @test qns+qn==qn+qns==QuantumNumbers('U',[qn1+qn,qn2+qn],[0,2,4],qnsindptr)
    @test qns-qn==QuantumNumbers('U',[qn1-qn,qn2-qn],[0,2,4],qnsindptr)
    @test qn-qns==QuantumNumbers('U',[qn-qn1,qn-qn2],[0,2,4],qnsindptr)

    qns1=QuantumNumbers('U',[+qn1,-qn2],[0,2,4],qnsindptr)
    qns2=QuantumNumbers('U',[-qn1,+qn2],[0,3,4],qnsindptr)
    @test ⊕((qns1,qns2),(+1,+1))==QuantumNumbers('G',[+qn1,-qn2,-qn1,+qn2],[2,2,3,1],qnscounts)
    @test ⊕((qns1,qns2),(+1,-1))==QuantumNumbers('G',[+qn1,-qn2,+qn1,-qn2],[2,2,3,1],qnscounts)
    @test ⊕((qns1,qns2),(-1,+1))==QuantumNumbers('G',[-qn1,+qn2,-qn1,+qn2],[2,2,3,1],qnscounts)
    @test ⊕((qns1,qns2),(-1,-1))==QuantumNumbers('G',[-qn1,+qn2,+qn1,-qn2],[2,2,3,1],qnscounts)
    @test ⊗((qns1,qns2),(+1,+1))==QuantumNumbers('G',[+qn1-qn1,+qn1+qn2,+qn1-qn1,+qn1+qn2,-qn2-qn1,-qn2+qn2,-qn2-qn1,-qn2+qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test ⊗((qns1,qns2),(+1,-1))==QuantumNumbers('G',[+qn1+qn1,+qn1-qn2,+qn1+qn1,+qn1-qn2,-qn2+qn1,-qn2-qn2,-qn2+qn1,-qn2-qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test ⊗((qns1,qns2),(-1,+1))==QuantumNumbers('G',[-qn1-qn1,-qn1+qn2,-qn1-qn1,-qn1+qn2,+qn2-qn1,+qn2+qn2,+qn2-qn1,+qn2+qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test ⊗((qns1,qns2),(-1,-1))==QuantumNumbers('G',[-qn1+qn1,-qn1-qn2,-qn1+qn1,-qn1-qn2,+qn2+qn1,+qn2-qn2,+qn2+qn1,+qn2-qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test kron((qns1,qns2),(+1,+1))==⊗((qns1,qns2),(+1,+1))
    @test kron((qns1,qns2),(+1,-1))==⊗((qns1,qns2),(+1,-1))
    @test kron((qns1,qns2),(-1,+1))==⊗((qns1,qns2),(-1,+1))
    @test kron((qns1,qns2),(-1,-1))==⊗((qns1,qns2),(-1,-1))
end

@testset "sort" begin
    qn1,qn2,qn3,qn4=CNZ4(3.0,2.0),CNZ4(4.0,1.0),CNZ4(4.0,2.0),CNZ4(3.0,3.0)
    oldqns=QuantumNumbers('G',[qn1,qn2,qn3,qn4,qn2,qn3,qn1,qn4],[2,1,3,4,1,1,2,3],qnscounts)
    newqns,permutation=sort(oldqns)
    @test newqns==QuantumNumbers('C',[qn1,qn4,qn2,qn3],[4,7,2,4],qnscounts)
    @test permutation==[1,2,13,14,7,8,9,10,15,16,17,3,11,4,5,6,12]
end

@testset "findall" begin
    qn1,qn2,qn3=CNZ4(2.0,3.0),CNZ4(1.0,2.0),CNZ4(0.0,0.0)
    qns=QuantumNumbers('C',[qn2,qn1],[2,3],qnscounts)
    @test findall(qns,qn1,qnscontents)==[2]
    @test findall(qns,qn2,qnscontents)==[1]
    @test findall(qns,qn3,qnscontents)==[]
    @test findall(qns,qn1,qnsexpansion)==[3,4,5]
    @test findall(qns,qn2,qnsexpansion)==[1,2]
    @test findall(qns,qn3,qnsexpansion)==[]

    qns=QuantumNumbers('G',[qn1,qn2,qn1],[2,2,1],qnscounts)
    @test findall(qns,qn1,qnscontents)==[1,3]
    @test findall(qns,qn2,qnscontents)==[2]
    @test findall(qns,qn3,qnscontents)==[]
    @test findall(qns,qn1,qnsexpansion)==[1,2,5]
    @test findall(qns,qn2,qnsexpansion)==[3,4]
    @test findall(qns,qn3,qnsexpansion)==[]
end

@testset "ukron" begin
    qn1,qn2=CNZ4(1.0,1.0),CNZ4(2.0,3.0)
    qns1=QuantumNumbers('U',[qn2,qn1],[4,5],qnscounts)
    qns2=QuantumNumbers('U',[qn1,qn2],[2,3],qnscounts)
    qns,records=ukron((qns1,qns2),(+1,-1))
    @test qns==QuantumNumbers('C',[qn1-qn2,CNZ4|>zero,qn2-qn1],[15,22,8],qnscounts)
    @test records==Dict(
        (qn1-qn2) => Dict((qn1,qn2)=>1:15),
        CNZ4|>zero => Dict((qn2,qn2)=>1:12,(qn1,qn1)=>13:22),
        (qn2-qn1) => Dict((qn2,qn1)=>1:8)
    )
end

@testset "expand" begin
    qn=CNZ4(3.0,2.0)
    qns=QuantumNumbers(qn,4)
    @test expand(qns,qnsindices)==[1,1,1,1]
    @test expand(qns,qnscontents)==[qn,qn,qn,qn]
end

@testset "decompose" begin
    qn1,qn2=CNZ4(0.0,0.0),CNZ4(1.0,1.0)
    qns=QuantumNumbers('U',[qn1,qn2],[0,1,2],qnsindptr)
    target=CNZ4(2.0,2.0)
    result=Set(((1,1,2,2),(1,2,1,2),(1,2,2,1),(2,1,1,2),(2,1,2,1),(2,2,1,1)))
    @test ⊆(Set(decompose((qns,qns,qns,qns),target,(1,1,1,1),qnsbruteforce,nmax=10)),result)
    @test ⊆(Set(decompose((qns,qns,qns,qns),target,(1,1,1,1),qnsmontecarlo,nmax=10)),result)
end

@testset "subset" begin
    qn1,qn2=CNZ4(1.0,2.0),CNZ4(3.0,0.0)
    qns=QuantumNumbers('G',[qn1,qn2,qn1,qn2],[1,2,3,4],qnscounts)
    @test subset(qns,qn1)==QuantumNumbers('G',[qn1,qn1],[1,3],qnscounts)
    @test subset(qns,qn2)==QuantumNumbers('G',[qn2,qn2],[2,4],qnscounts)
end

@testset "reorder" begin
    qn1,qn2,qn3=CNZ4(1.0,2.0),CNZ4(3.0,0.0),CNZ4(4.0,1.0)
    qns=QuantumNumbers('G',[qn1,qn2,qn3],[2,3,4],qnscounts)
    @test reorder(qns,[3,2,1],qnscontents)==QuantumNumbers('G',[qn3,qn2,qn1],[4,3,2],qnscounts)
    @test reorder(qns,[4,6,9,8],qnsexpansion)==QuantumNumbers('G',[qn2,qn3,qn3,qn3],[1,1,1,1],qnscounts)
end

@testset "toordereddict" begin
    qn1,qn2=CNZ4(1.0,2.0),CNZ4(2.0,3.0)
    qns=QuantumNumbers('U',[qn1,qn2],[2,3],qnscounts)
    @test toordereddict(qns,qnsindptr)==OrderedDict(qn1=>1:2,qn2=>3:5)
    @test toordereddict(qns,qnscounts)==OrderedDict(qn1=>2,qn2=>3)
end
