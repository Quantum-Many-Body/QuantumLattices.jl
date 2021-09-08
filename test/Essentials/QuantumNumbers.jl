using Test
using Base.Iterators: Iterators
using DataStructures: OrderedDict
using QuantumLattices.Essentials.QuantumNumbers
using QuantumLattices.Interfaces: ⊕, ⊗, dimension, expand, permute, decompose
using QuantumLattices.Prerequisites.Traits: contentnames, getcontent

import QuantumLattices.Essentials.QuantumNumbers: periods

@abeliannumber "CN" (:N,) (Inf,)
@abeliannumber "Z4" (:Z,) (4,)
@abeliannumber "CNZ4" (:N, :Z) (Inf, 4)

@testset "regularize" begin
    @test regularize(CNZ4, [1.5, 5.0]) == [1.5, 1.0]
    @test regularize(CNZ4, [4.0, -1.0]) == [4.0, 3.0]
    @test regularize(CNZ4, [1.5 4.0; 5.0 -1.0]) == [1.5 4.0; 1.0 3.0]
end

@testset "AbelianNumbers" begin
    @test contentnames(AbelianNumbers) == (:form, :table, :indptr)

    qn1, qn2 = CNZ4(+1.0, +3.0), CNZ4(-1.0, +1.0)
    qns1 = AbelianNumbers('U', [qn1, qn2], [0, 2, 5], :indptr)
    qns2 = AbelianNumbers('U', [qn1, qn2], [2, 3], :counts)
    @test isequal(qns1, qns2)
    qns3 = AbelianNumbers(OrderedDict(qn1=>2, qn2=>3))
    qns4 = AbelianNumbers(OrderedDict(qn1=>1:2, qn2=>3:5))
    @test qns1 == qns2 == qns3 == qns4

    qns = AbelianNumbers('U', [qn1, qn2], [0, 3, 5], :indptr)
    @test getcontent(qns, Val(:table)) == qns.contents
    @test qns|>dimension == 5
    @test qns|>string == "QNS(2, 5)"
    @test qns|>size == (2,)
    @test qns|>eltype == CNZ4
    @test qns|>typeof|>eltype == CNZ4
    @test qns|>issorted == false
    @test repr(qns) == "QNS(CNZ4(1.0, 3.0)=>1:3, CNZ4(-1.0, 1.0)=>4:5)"
    @test qns[1] == qn1
    @test qns[2] == qn2
    @test qns[1:2] == AbelianNumbers('U', [qn1, qn2], [0, 3, 5], :indptr)
    @test qns[[2, 1]] == AbelianNumbers('G', [qn2, qn1], [0, 2, 5], :indptr)
    @test qns|>collect == [qn1, qn2]
    @test qns|>Iterators.reverse|>collect == [qn2, qn1]
    @test qns|>keys|>collect == [qn1, qn2]
    @test values(qns, :indptr)|>collect == [1:3, 4:5]
    @test values(qns, :counts)|>collect == [3, 2]
    @test pairs(qns, :indptr)|>collect == [qn1=>1:3, qn2=>4:5]
    @test pairs(qns, :counts)|>collect == [qn1=>3, qn2=>2]
end

@testset "OrderedDict" begin
    qn1, qn2 = CNZ4(1.0, 2.0), CNZ4(2.0, 3.0)
    qns = AbelianNumbers('U', [qn1, qn2], [2, 3], :counts)
    @test OrderedDict(qns, :indptr) == OrderedDict(qn1=>1:2, qn2=>3:5)
    @test OrderedDict(qns, :counts) == OrderedDict(qn1=>2, qn2=>3)
end

@testset "arithmetic" begin
    qn = CNZ4(1.0, 3.0)
    @test qn|>dimension == 1
    @test qn|>typeof|>dimension == 1
    @test +qn == CNZ4(+1.0, +3.0)
    @test -qn == CNZ4(-1.0, +1.0)
    @test qn*4 == 4*qn == CNZ4(4.0, 0.0)
    @test qn^3 == CNZ4(3.0, 1.0)

    qn1, qn2 = CNZ4(1.0, 2.0), CNZ4(2.0, 3.0)
    @test qn1+qn2 == CNZ4(3.0, 1.0)
    @test qn1-qn2 == CNZ4(-1.0, 3.0)
    @test kron(qn1, qn2, signs=(+1, +1)) == CNZ4(+3.0, 1.0)
    @test kron(qn1, qn2, signs=(+1, -1)) == CNZ4(-1.0, 3.0)
    @test kron(qn1, qn2, signs=(-1, +1)) == CNZ4(+1.0, 1.0)
    @test kron(qn1, qn2, signs=(-1, -1)) == CNZ4(-3.0, 3.0)
    @test union(qn1, qn2, signs=(+1, +1)) == AbelianNumbers('G', [+qn1, +qn2], [0, 1, 2], :indptr)
    @test union(qn1, qn2, signs=(+1, -1)) == AbelianNumbers('G', [+qn1, -qn2], [0, 1, 2], :indptr)
    @test union(qn1, qn2, signs=(-1, +1)) == AbelianNumbers('G', [-qn1, +qn2], [0, 1, 2], :indptr)
    @test union(qn1, qn2, signs=(-1, -1)) == AbelianNumbers('G', [-qn1, -qn2], [0, 1, 2], :indptr)
    @test ⊗(qn1, qn2) == kron(qn1, qn2)
    @test ⊕(qn1, qn2) == union(qn1, qn2)

    qns = AbelianNumbers('U', [qn1, qn2], [0, 2, 4], :indptr)
    @test +qns == AbelianNumbers('U', [+qn1, +qn2], [0, 2, 4], :indptr)
    @test -qns == AbelianNumbers('U', [-qn1, -qn2], [0, 2, 4], :indptr)
    @test qns*3 == 3*qns == AbelianNumbers('U', [qn1*3, qn2*3], [0, 2, 4], :indptr)
    @test qns^2 == AbelianNumbers('G', [qn1+qn1, qn1+qn2, qn1+qn1, qn1+qn2, qn2+qn1, qn2+qn2, qn2+qn1, qn2+qn2], [2, 2, 2, 2, 2, 2, 2, 2], :counts)

    @test qns+qn == qn+qns == AbelianNumbers('U', [qn1+qn, qn2+qn], [0, 2, 4], :indptr)
    @test qns-qn == AbelianNumbers('U', [qn1-qn, qn2-qn], [0, 2, 4], :indptr)
    @test qn-qns == AbelianNumbers('U', [qn-qn1, qn-qn2], [0, 2, 4], :indptr)

    qns1 = AbelianNumbers('U', [+qn1, -qn2], [0, 2, 4], :indptr)
    qns2 = AbelianNumbers('U', [-qn1, +qn2], [0, 3, 4], :indptr)
    @test union(qns1, qns2, signs=(+1, +1)) == AbelianNumbers('G', [+qn1, -qn2, -qn1, +qn2], [2, 2, 3, 1], :counts)
    @test union(qns1, qns2, signs=(+1, -1)) == AbelianNumbers('G', [+qn1, -qn2, +qn1, -qn2], [2, 2, 3, 1], :counts)
    @test union(qns1, qns2, signs=(-1, +1)) == AbelianNumbers('G', [-qn1, +qn2, -qn1, +qn2], [2, 2, 3, 1], :counts)
    @test union(qns1, qns2, signs=(-1, -1)) == AbelianNumbers('G', [-qn1, +qn2, +qn1, -qn2], [2, 2, 3, 1], :counts)
    @test kron(qns1, qns2, signs=(+1, +1)) == AbelianNumbers('G', [+qn1-qn1, +qn1+qn2, +qn1-qn1, +qn1+qn2, -qn2-qn1, -qn2+qn2, -qn2-qn1, -qn2+qn2], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns1, qns2, signs=(+1, -1)) == AbelianNumbers('G', [+qn1+qn1, +qn1-qn2, +qn1+qn1, +qn1-qn2, -qn2+qn1, -qn2-qn2, -qn2+qn1, -qn2-qn2], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns1, qns2, signs=(-1, +1)) == AbelianNumbers('G', [-qn1-qn1, -qn1+qn2, -qn1-qn1, -qn1+qn2, +qn2-qn1, +qn2+qn2, +qn2-qn1, +qn2+qn2], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns1, qns2, signs=(-1, -1)) == AbelianNumbers('G', [-qn1+qn1, -qn1-qn2, -qn1+qn1, -qn1-qn2, +qn2+qn1, +qn2-qn2, +qn2+qn1, +qn2-qn2], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test ⊕(qns1, qns2) == union(qns1, qns2)
    @test ⊗(qns1, qns2) == kron(qns1, qns2)
end

@testset "sort" begin
    qn1, qn2, qn3, qn4 = CNZ4(3.0, 2.0), CNZ4(4.0, 1.0), CNZ4(4.0, 2.0), CNZ4(3.0, 3.0)
    oldqns = AbelianNumbers('G', [qn1, qn2, qn3, qn4, qn2, qn3, qn1, qn4], [2, 1, 3, 4, 1, 1, 2, 3], :counts)
    newqns, permutation = sort(oldqns)
    @test newqns == AbelianNumbers('C', [qn1, qn4, qn2, qn3], [4, 7, 2, 4], :counts)
    @test permutation == [1, 2, 13, 14, 7, 8, 9, 10, 15, 16, 17, 3, 11, 4, 5, 6, 12]
end

@testset "findall" begin
    qn1, qn2, qn3 = CNZ4(2.0, 3.0), CNZ4(1.0, 2.0), CNZ4(0.0, 0.0)
    qns = AbelianNumbers('C', [qn2, qn1], [2, 3], :counts)
    @test findall(qn1, qns, :compression) == [2]
    @test findall(qn2, qns, :compression) == [1]
    @test findall(qn3, qns, :compression) == []
    @test findall(qn1, qns, :expansion) == [3, 4, 5]
    @test findall(qn2, qns, :expansion) == [1, 2]
    @test findall(qn3, qns, :expansion) == []

    qns = AbelianNumbers('G', [qn1, qn2, qn1], [2, 2, 1], :counts)
    @test findall(qn1, qns, :compression) == [1, 3]
    @test findall(qn2, qns, :compression) == [2]
    @test findall(qn3, qns, :compression) == []
    @test findall(qn1, qns, :expansion) == [1, 2, 5]
    @test findall(qn2, qns, :expansion) == [3, 4]
    @test findall(qn3, qns, :expansion) == []
end

@testset "filter" begin
    qn1, qn2 = CNZ4(1.0, 2.0), CNZ4(3.0, 0.0)
    qns = AbelianNumbers('G', [qn1, qn2, qn1, qn2], [1, 2, 3, 4], :counts)
    @test filter((qn1, qn2), qns) == AbelianNumbers('G', [qn1, qn2, qn1, qn2], [1, 2, 3, 4], :counts)
    @test filter(qn2, qns) == AbelianNumbers('G', [qn2, qn2], [2, 4], :counts)
end

@testset "prod" begin
    qn1, qn2 = CNZ4(1.0, 1.0), CNZ4(2.0, 3.0)
    qns1 = AbelianNumbers('U', [qn2, qn1], [4, 5], :counts)
    qns2 = AbelianNumbers('U', [qn1, qn2], [2, 3], :counts)
    qns, records = prod(qns1, qns2, signs=(+1, -1))
    @test qns == AbelianNumbers('C', [qn1-qn2, zero(CNZ4), qn2-qn1], [15, 22, 8], :counts)
    @test records == Dict(  (qn1-qn2) =>  Dict((qn1, qn2)=>1:15),
                            zero(CNZ4) => Dict((qn2, qn2)=>1:12, (qn1, qn1)=>13:22),
                            (qn2-qn1) =>  Dict((qn2, qn1)=>1:8)
                            )
end

@testset "expand" begin
    qn = CNZ4(3.0, 2.0)
    qns = AbelianNumbers(qn, 4)
    @test expand(qns, :indices) == [1, 1, 1, 1]
    @test expand(qns, :contents) == [qn, qn, qn, qn]
end

@testset "decompose" begin
    qn1, qn2 = CNZ4(0.0, 0.0), CNZ4(1.0, 1.0)
    qns = AbelianNumbers('U', [qn1, qn2], [0, 1, 2], :indptr)
    target = CNZ4(2.0, 2.0)
    result = Set(((1, 1, 2, 2), (1, 2, 1, 2), (1, 2, 2, 1), (2, 1, 1, 2), (2, 1, 2, 1), (2, 2, 1, 1)))
    @test ⊆(Set(decompose(target, qns, qns, qns, qns; signs=(1, 1, 1, 1), nmax=10, method=:bruteforce)), result)
    @test ⊆(Set(decompose(target, qns, qns, qns, qns; signs=(1, 1, 1, 1), nmax=10, method=:montecarlo)), result)
end

@testset "permute" begin
    qn1, qn2, qn3 = CNZ4(1.0, 2.0), CNZ4(3.0, 0.0), CNZ4(4.0, 1.0)
    qns = AbelianNumbers('G', [qn1, qn2, qn3], [2, 3, 4], :counts)
    @test permute(qns, [3, 2, 1], :compression) == AbelianNumbers('G', [qn3, qn2, qn1], [4, 3, 2], :counts)
    @test permute(qns, [4, 6, 9, 8], :expansion) == AbelianNumbers('G', [qn2, qn3, qn3, qn3], [1, 1, 1, 1], :counts)
end

@testset "ConcreteAbelianNumbers" begin
    @test SQN|>fieldnames == (:Sz,)
    @test SQN|>periods == (Inf,)
    @test PQN|>fieldnames == (:N,)
    @test PQN|>periods == (Inf,)
    @test SPQN|>fieldnames == (:N, :Sz)
    @test SPQN|>periods == (Inf, Inf)
    @test SQNS(0.5) == AbelianNumbers('C', [SQN(-0.5), SQN(0.5)], [0, 1, 2], :indptr)
    @test PQNS(1.0) == AbelianNumbers('C', [PQN(0.0), PQN(1.0)], [0, 1, 2], :indptr)
    @test SPQNS(0.5) == AbelianNumbers('C', [SPQN(0.0, 0.0), SPQN(1.0, -0.5), SPQN(1.0, 0.5), SPQN(2.0, 0.0)], [0, 1, 2, 3, 4], :indptr)

    @test Momentum1D{10} |> periods == (10,)
    @test Momentum2D{10, 15} |> periods == (10, 15)
    @test Momentum3D{10, 15, 20} |> periods == (10, 15, 20)

    @test Momentum1D{10}(1) == Momentum1D{10}(11) == Momentum1D{10}(-9)
    @test Momentum2D{10}(1, 1) == Momentum2D{10}(11, 11) == Momentum2D{10}(-9, -9)
    @test Momentum3D{10}(1, 1, 1) == Momentum3D{10}(11, 11, 11) == Momentum3D{10}(-9, -9, -9)
    @test Momentum2D{10, 20}(1, 1) == Momentum2D{10, 20}(11, 21) == Momentum2D{10, 20}(-9, -19)
    @test Momentum3D{10, 20, 30}(1, 1, 1) == Momentum3D{10, 20, 30}(11, 21, 31) == Momentum3D{10, 20, 30}(-9, -19, -29)
end