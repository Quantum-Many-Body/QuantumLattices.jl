using Test
using DataStructures: OrderedDict
using Printf: @sprintf
using QuantumLattices.Mathematics.QuantumNumbers
using QuantumLattices.Interfaces: ⊕,⊗,dimension,expand,permute,decompose,regularize
import Base.Iterators: Iterators

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
    @test SQNS(0.5)==AbelianNumbers('C',[SQN(-0.5),SQN(0.5)],[0,1,2],qnsindptr)
end

@testset "PQNS" begin
    @test PQNS(1.0)==AbelianNumbers('C',[PQN(0.0),PQN(1.0)],[0,1,2],qnsindptr)
end

@testset "SzPQNS" begin
    @test SzPQNS(0.5)==AbelianNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,0.5)],[0,1,2],qnsindptr)
end

@testset "SPQNS" begin
    @test SPQNS(0.5)==AbelianNumbers('C',[SPQN(0.0,0.0),SPQN(1.0,-0.5),SPQN(1.0,0.5),SPQN(2.0,0.0)],[0,1,2,3,4],qnsindptr)
end

@testset "Z2QNS" begin
    @test Z2QNS()==AbelianNumbers('C',[Z2QN(0.0),Z2QN(1.0)],[0,1,2],qnsindptr)
end

@abeliannumber "CN" (:N,) (Inf,)
@abeliannumber "Z4" (:Z,) (4,)
@abeliannumber "CNZ4" (:N,:Z) (Inf,4)

@testset "regularize" begin
    @test regularize(CNZ4,[1.5,5.0])==[1.5,1.0]
    @test regularize(CNZ4,[4.0,-1.0])==[4.0,3.0]
    @test regularize(CNZ4,[1.5 4.0; 5.0 -1.0])==[1.5 4.0; 1.0 3.0]
end

@testset "AbelianNumbers" begin
    qn1,qn2=CNZ4(+1.0,+3.0),CNZ4(-1.0,+1.0)
    qns1=AbelianNumbers('U',[qn1,qn2],[0,2,5],qnsindptr)
    qns2=AbelianNumbers('U',[qn1,qn2],[2,3],qnscounts)
    @test isequal(qns1,qns2)
    qns3=AbelianNumbers(OrderedDict(qn1=>2,qn2=>3))
    qns4=AbelianNumbers(OrderedDict(qn1=>1:2,qn2=>3:5))
    @test qns1==qns2==qns3==qns4

    qns=AbelianNumbers('U',[qn1,qn2],[0,3,5],qnsindptr)
    @test qns|>dimension==5
    @test qns|>string=="QNS(2,5)"
    @test qns|>length==2
    @test qns|>eltype==CNZ4
    @test qns|>typeof|>eltype==CNZ4
    @test @sprintf("%s",qns)=="QNS(CNZ4(1.0,3.0)=>1:3,CNZ4(-1.0,1.0)=>4:5)"
    @test qns[1]==qn1
    @test qns[2]==qn2
    @test qns[1:2]==AbelianNumbers('G',[qn1,qn2],[0,3,5],qnsindptr)
    @test qns[[2,1]]==AbelianNumbers('G',[qn2,qn1],[0,2,5],qnsindptr)
    @test qns|>collect==[qn1,qn2]
    @test qns|>Iterators.reverse|>collect==[qn2,qn1]
    @test qns|>keys|>collect==[qn1,qn2]
    @test values(qns,qnsindptr)|>collect==[1:3,4:5]
    @test values(qns,qnscounts)|>collect==[3,2]
    @test pairs(qns,qnsindptr)|>collect==[qn1=>1:3,qn2=>4:5]
    @test pairs(qns,qnscounts)|>collect==[qn1=>3,qn2=>2]
end

@testset "toordereddict" begin
    qn1,qn2=CNZ4(1.0,2.0),CNZ4(2.0,3.0)
    qns=AbelianNumbers('U',[qn1,qn2],[2,3],qnscounts)
    @test toordereddict(qns,qnsindptr)==OrderedDict(qn1=>1:2,qn2=>3:5)
    @test toordereddict(qns,qnscounts)==OrderedDict(qn1=>2,qn2=>3)
end

@testset "arithmetic" begin
    qn=CNZ4(1.0,3.0)
    @test qn|>dimension==1
    @test qn|>typeof|>dimension==1
    @test +qn==CNZ4(+1.0,+3.0)
    @test -qn==CNZ4(-1.0,+1.0)
    @test qn*4==4*qn==CNZ4(4.0,0.0)
    @test qn^3==CNZ4(3.0,1.0)
    @test kron(CNZ4,CN(1.0),Z4(3.0))==qn

    qn1,qn2=CNZ4(1.0,2.0),CNZ4(2.0,3.0)
    @test qn1+qn2==CNZ4(3.0,1.0)
    @test qn1-qn2==CNZ4(-1.0,3.0)
    @test kron(qn1,qn2,signs=(+1,+1))==CNZ4(+3.0,1.0)
    @test kron(qn1,qn2,signs=(+1,-1))==CNZ4(-1.0,3.0)
    @test kron(qn1,qn2,signs=(-1,+1))==CNZ4(+1.0,1.0)
    @test kron(qn1,qn2,signs=(-1,-1))==CNZ4(-3.0,3.0)
    @test union(qn1,qn2,signs=(+1,+1))==AbelianNumbers('G',[+qn1,+qn2],[0,1,2],qnsindptr)
    @test union(qn1,qn2,signs=(+1,-1))==AbelianNumbers('G',[+qn1,-qn2],[0,1,2],qnsindptr)
    @test union(qn1,qn2,signs=(-1,+1))==AbelianNumbers('G',[-qn1,+qn2],[0,1,2],qnsindptr)
    @test union(qn1,qn2,signs=(-1,-1))==AbelianNumbers('G',[-qn1,-qn2],[0,1,2],qnsindptr)
    @test ⊗(qn1,qn2)==kron(qn1,qn2)
    @test ⊕(qn1,qn2)==union(qn1,qn2)

    qns=AbelianNumbers('U',[qn1,qn2],[0,2,4],qnsindptr)
    @test +qns==AbelianNumbers('U',[+qn1,+qn2],[0,2,4],qnsindptr)
    @test -qns==AbelianNumbers('U',[-qn1,-qn2],[0,2,4],qnsindptr)
    @test qns*3==3*qns==AbelianNumbers('U',[qn1*3,qn2*3],[0,2,4],qnsindptr)
    @test qns^2==AbelianNumbers('G',[qn1+qn1,qn1+qn2,qn1+qn1,qn1+qn2,qn2+qn1,qn2+qn2,qn2+qn1,qn2+qn2],[2,2,2,2,2,2,2,2],qnscounts)

    @test qns+qn==qn+qns==AbelianNumbers('U',[qn1+qn,qn2+qn],[0,2,4],qnsindptr)
    @test qns-qn==AbelianNumbers('U',[qn1-qn,qn2-qn],[0,2,4],qnsindptr)
    @test qn-qns==AbelianNumbers('U',[qn-qn1,qn-qn2],[0,2,4],qnsindptr)

    qns1=AbelianNumbers('U',[+qn1,-qn2],[0,2,4],qnsindptr)
    qns2=AbelianNumbers('U',[-qn1,+qn2],[0,3,4],qnsindptr)
    @test union(qns1,qns2,signs=(+1,+1))==AbelianNumbers('G',[+qn1,-qn2,-qn1,+qn2],[2,2,3,1],qnscounts)
    @test union(qns1,qns2,signs=(+1,-1))==AbelianNumbers('G',[+qn1,-qn2,+qn1,-qn2],[2,2,3,1],qnscounts)
    @test union(qns1,qns2,signs=(-1,+1))==AbelianNumbers('G',[-qn1,+qn2,-qn1,+qn2],[2,2,3,1],qnscounts)
    @test union(qns1,qns2,signs=(-1,-1))==AbelianNumbers('G',[-qn1,+qn2,+qn1,-qn2],[2,2,3,1],qnscounts)
    @test kron(qns1,qns2,signs=(+1,+1))==AbelianNumbers('G',[+qn1-qn1,+qn1+qn2,+qn1-qn1,+qn1+qn2,-qn2-qn1,-qn2+qn2,-qn2-qn1,-qn2+qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test kron(qns1,qns2,signs=(+1,-1))==AbelianNumbers('G',[+qn1+qn1,+qn1-qn2,+qn1+qn1,+qn1-qn2,-qn2+qn1,-qn2-qn2,-qn2+qn1,-qn2-qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test kron(qns1,qns2,signs=(-1,+1))==AbelianNumbers('G',[-qn1-qn1,-qn1+qn2,-qn1-qn1,-qn1+qn2,+qn2-qn1,+qn2+qn2,+qn2-qn1,+qn2+qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test kron(qns1,qns2,signs=(-1,-1))==AbelianNumbers('G',[-qn1+qn1,-qn1-qn2,-qn1+qn1,-qn1-qn2,+qn2+qn1,+qn2-qn2,+qn2+qn1,+qn2-qn2],[3,1,3,1,3,1,3,1],qnscounts)
    @test ⊕(qns1,qns2)==union(qns1,qns2)
    @test ⊗(qns1,qns2)==kron(qns1,qns2)
end

@testset "sort" begin
    qn1,qn2,qn3,qn4=CNZ4(3.0,2.0),CNZ4(4.0,1.0),CNZ4(4.0,2.0),CNZ4(3.0,3.0)
    oldqns=AbelianNumbers('G',[qn1,qn2,qn3,qn4,qn2,qn3,qn1,qn4],[2,1,3,4,1,1,2,3],qnscounts)
    newqns,permutation=sort(oldqns)
    @test newqns==AbelianNumbers('C',[qn1,qn4,qn2,qn3],[4,7,2,4],qnscounts)
    @test permutation==[1,2,13,14,7,8,9,10,15,16,17,3,11,4,5,6,12]
end

@testset "findall" begin
    qn1,qn2,qn3=CNZ4(2.0,3.0),CNZ4(1.0,2.0),CNZ4(0.0,0.0)
    qns=AbelianNumbers('C',[qn2,qn1],[2,3],qnscounts)
    @test findall(qn1,qns,qnscompression)==[2]
    @test findall(qn2,qns,qnscompression)==[1]
    @test findall(qn3,qns,qnscompression)==[]
    @test findall(qn1,qns,qnsexpansion)==[3,4,5]
    @test findall(qn2,qns,qnsexpansion)==[1,2]
    @test findall(qn3,qns,qnsexpansion)==[]

    qns=AbelianNumbers('G',[qn1,qn2,qn1],[2,2,1],qnscounts)
    @test findall(qn1,qns,qnscompression)==[1,3]
    @test findall(qn2,qns,qnscompression)==[2]
    @test findall(qn3,qns,qnscompression)==[]
    @test findall(qn1,qns,qnsexpansion)==[1,2,5]
    @test findall(qn2,qns,qnsexpansion)==[3,4]
    @test findall(qn3,qns,qnsexpansion)==[]
end

@testset "filter" begin
    qn1,qn2=CNZ4(1.0,2.0),CNZ4(3.0,0.0)
    qns=AbelianNumbers('G',[qn1,qn2,qn1,qn2],[1,2,3,4],qnscounts)
    @test filter((qn1,qn2),qns)==AbelianNumbers('G',[qn1,qn2,qn1,qn2],[1,2,3,4],qnscounts)
    @test filter(qn2,qns)==AbelianNumbers('G',[qn2,qn2],[2,4],qnscounts)
end

@testset "ukron" begin
    qn1,qn2=CNZ4(1.0,1.0),CNZ4(2.0,3.0)
    qns1=AbelianNumbers('U',[qn2,qn1],[4,5],qnscounts)
    qns2=AbelianNumbers('U',[qn1,qn2],[2,3],qnscounts)
    qns,records=ukron(qns1,qns2,signs=(+1,-1))
    @test qns==AbelianNumbers('C',[qn1-qn2,CNZ4|>zero,qn2-qn1],[15,22,8],qnscounts)
    @test records==Dict(
        (qn1-qn2) => Dict((qn1,qn2)=>1:15),
        CNZ4|>zero => Dict((qn2,qn2)=>1:12,(qn1,qn1)=>13:22),
        (qn2-qn1) => Dict((qn2,qn1)=>1:8)
    )
end

@testset "expand" begin
    qn=CNZ4(3.0,2.0)
    qns=AbelianNumbers(qn,4)
    @test expand(qns,qnsindices)==[1,1,1,1]
    @test expand(qns,qnscontents)==[qn,qn,qn,qn]
end

@testset "decompose" begin
    qn1,qn2=CNZ4(0.0,0.0),CNZ4(1.0,1.0)
    qns=AbelianNumbers('U',[qn1,qn2],[0,1,2],qnsindptr)
    target=CNZ4(2.0,2.0)
    result=Set(((1,1,2,2),(1,2,1,2),(1,2,2,1),(2,1,1,2),(2,1,2,1),(2,2,1,1)))
    @test ⊆(Set(decompose((qns,qns,qns,qns),target,(1,1,1,1),qnsbruteforce,nmax=10)),result)
    @test ⊆(Set(decompose((qns,qns,qns,qns),target,(1,1,1,1),qnsmontecarlo,nmax=10)),result)
end

@testset "permute" begin
    qn1,qn2,qn3=CNZ4(1.0,2.0),CNZ4(3.0,0.0),CNZ4(4.0,1.0)
    qns=AbelianNumbers('G',[qn1,qn2,qn3],[2,3,4],qnscounts)
    @test permute(qns,[3,2,1],qnscompression)==AbelianNumbers('G',[qn3,qn2,qn1],[4,3,2],qnscounts)
    @test permute(qns,[4,6,9,8],qnsexpansion)==AbelianNumbers('G',[qn2,qn3,qn3,qn3],[1,1,1,1],qnscounts)
end
