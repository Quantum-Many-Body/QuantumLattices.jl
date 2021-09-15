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

    qn₁, qn₂ = CNZ4(+1.0, +3.0), CNZ4(-1.0, +1.0)
    qns₁ = AbelianNumbers('U', [qn₁, qn₂], [0, 2, 5], :indptr)
    qns₂ = AbelianNumbers('U', [qn₁, qn₂], [2, 3], :counts)
    @test isequal(qns₁, qns₂)
    qns₃ = AbelianNumbers(OrderedDict(qn₁=>2, qn₂=>3))
    qns₄ = AbelianNumbers(OrderedDict(qn₁=>1:2, qn₂=>3:5))
    @test qns₁ == qns₂ == qns₃ == qns₄

    qns = AbelianNumbers('U', [qn₁, qn₂])
    @test qns == AbelianNumbers('U', [qn₁, qn₂], [1, 1], :counts)

    qns = AbelianNumbers('U', [qn₁, qn₂], [0, 3, 5], :indptr)
    @test getcontent(qns, Val(:table)) == qns.contents
    @test qns|>dimension == 5
    @test qns|>string == "QNS(2, 5)"
    @test qns|>size == (2,)
    @test qns|>eltype == CNZ4
    @test qns|>typeof|>eltype == CNZ4
    @test qns|>issorted == false
    @test repr(qns) == "QNS(CNZ4(1.0, 3.0)=>1:3, CNZ4(-1.0, 1.0)=>4:5)"
    @test qns[1] == qn₁
    @test qns[2] == qn₂
    @test qns[1:2] == AbelianNumbers('U', [qn₁, qn₂], [0, 3, 5], :indptr)
    @test qns[[2, 1]] == AbelianNumbers('G', [qn₂, qn₁], [0, 2, 5], :indptr)
    @test qns|>collect == [qn₁, qn₂]
    @test qns|>Iterators.reverse|>collect == [qn₂, qn₁]
    @test qns|>keys|>collect == [qn₁, qn₂]
    @test values(qns, :indptr)|>collect == [1:3, 4:5]
    @test values(qns, :counts)|>collect == [3, 2]
    @test pairs(qns, :indptr)|>collect == [qn₁=>1:3, qn₂=>4:5]
    @test pairs(qns, :counts)|>collect == [qn₁=>3, qn₂=>2]
end

@testset "OrderedDict" begin
    qn₁, qn₂ = CNZ4(1.0, 2.0), CNZ4(2.0, 3.0)
    qns = AbelianNumbers('U', [qn₁, qn₂], [2, 3], :counts)
    @test OrderedDict(qns, :indptr) == OrderedDict(qn₁=>1:2, qn₂=>3:5)
    @test OrderedDict(qns, :counts) == OrderedDict(qn₁=>2, qn₂=>3)
end

@testset "arithmetic" begin
    qn = CNZ4(1.0, 3.0)
    @test qn|>dimension == 1
    @test qn|>typeof|>dimension == 1
    @test +qn == CNZ4(+1.0, +3.0)
    @test -qn == CNZ4(-1.0, +1.0)
    @test qn*4 == 4*qn == CNZ4(4.0, 0.0)
    @test qn^3 == CNZ4(3.0, 1.0)

    qn₁, qn₂ = CNZ4(1.0, 2.0), CNZ4(2.0, 3.0)
    @test qn₁+qn₂ == CNZ4(3.0, 1.0)
    @test qn₁-qn₂ == CNZ4(-1.0, 3.0)
    @test kron(qn₁, qn₂, signs=(+1, +1)) == CNZ4(+3.0, 1.0)
    @test kron(qn₁, qn₂, signs=(+1, -1)) == CNZ4(-1.0, 3.0)
    @test kron(qn₁, qn₂, signs=(-1, +1)) == CNZ4(+1.0, 1.0)
    @test kron(qn₁, qn₂, signs=(-1, -1)) == CNZ4(-3.0, 3.0)
    @test union(qn₁, qn₂, signs=(+1, +1)) == AbelianNumbers('G', [+qn₁, +qn₂], [0, 1, 2], :indptr)
    @test union(qn₁, qn₂, signs=(+1, -1)) == AbelianNumbers('G', [+qn₁, -qn₂], [0, 1, 2], :indptr)
    @test union(qn₁, qn₂, signs=(-1, +1)) == AbelianNumbers('G', [-qn₁, +qn₂], [0, 1, 2], :indptr)
    @test union(qn₁, qn₂, signs=(-1, -1)) == AbelianNumbers('G', [-qn₁, -qn₂], [0, 1, 2], :indptr)
    @test ⊗(qn₁, qn₂) == kron(qn₁, qn₂)
    @test ⊕(qn₁, qn₂) == union(qn₁, qn₂)

    qns = AbelianNumbers('U', [qn₁, qn₂], [0, 2, 4], :indptr)
    @test +qns == AbelianNumbers('U', [+qn₁, +qn₂], [0, 2, 4], :indptr)
    @test -qns == AbelianNumbers('U', [-qn₁, -qn₂], [0, 2, 4], :indptr)
    @test qns*3 == 3*qns == AbelianNumbers('U', [qn₁*3, qn₂*3], [0, 2, 4], :indptr)
    @test qns^2 == AbelianNumbers('G', [qn₁+qn₁, qn₁+qn₂, qn₁+qn₁, qn₁+qn₂, qn₂+qn₁, qn₂+qn₂, qn₂+qn₁, qn₂+qn₂], [2, 2, 2, 2, 2, 2, 2, 2], :counts)

    @test qns+qn == qn+qns == AbelianNumbers('U', [qn₁+qn, qn₂+qn], [0, 2, 4], :indptr)
    @test qns-qn == AbelianNumbers('U', [qn₁-qn, qn₂-qn], [0, 2, 4], :indptr)
    @test qn-qns == AbelianNumbers('U', [qn-qn₁, qn-qn₂], [0, 2, 4], :indptr)

    qns₁ = AbelianNumbers('U', [+qn₁, -qn₂], [0, 2, 4], :indptr)
    qns₂ = AbelianNumbers('U', [-qn₁, +qn₂], [0, 3, 4], :indptr)
    @test union(qns₁, qns₂, signs=(+1, +1)) == AbelianNumbers('G', [+qn₁, -qn₂, -qn₁, +qn₂], [2, 2, 3, 1], :counts)
    @test union(qns₁, qns₂, signs=(+1, -1)) == AbelianNumbers('G', [+qn₁, -qn₂, +qn₁, -qn₂], [2, 2, 3, 1], :counts)
    @test union(qns₁, qns₂, signs=(-1, +1)) == AbelianNumbers('G', [-qn₁, +qn₂, -qn₁, +qn₂], [2, 2, 3, 1], :counts)
    @test union(qns₁, qns₂, signs=(-1, -1)) == AbelianNumbers('G', [-qn₁, +qn₂, +qn₁, -qn₂], [2, 2, 3, 1], :counts)
    @test kron(qns₁, qns₂, signs=(+1, +1)) == AbelianNumbers('G', [+qn₁-qn₁, +qn₁+qn₂, +qn₁-qn₁, +qn₁+qn₂, -qn₂-qn₁, -qn₂+qn₂, -qn₂-qn₁, -qn₂+qn₂], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns₁, qns₂, signs=(+1, -1)) == AbelianNumbers('G', [+qn₁+qn₁, +qn₁-qn₂, +qn₁+qn₁, +qn₁-qn₂, -qn₂+qn₁, -qn₂-qn₂, -qn₂+qn₁, -qn₂-qn₂], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns₁, qns₂, signs=(-1, +1)) == AbelianNumbers('G', [-qn₁-qn₁, -qn₁+qn₂, -qn₁-qn₁, -qn₁+qn₂, +qn₂-qn₁, +qn₂+qn₂, +qn₂-qn₁, +qn₂+qn₂], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test kron(qns₁, qns₂, signs=(-1, -1)) == AbelianNumbers('G', [-qn₁+qn₁, -qn₁-qn₂, -qn₁+qn₁, -qn₁-qn₂, +qn₂+qn₁, +qn₂-qn₂, +qn₂+qn₁, +qn₂-qn₂], [3, 1, 3, 1, 3, 1, 3, 1], :counts)
    @test ⊕(qns₁, qns₂) == union(qns₁, qns₂)
    @test ⊗(qns₁, qns₂) == kron(qns₁, qns₂)
end

@testset "sort" begin
    qn₁, qn₂, qn₃, qn₄ = CNZ4(3.0, 2.0), CNZ4(4.0, 1.0), CNZ4(4.0, 2.0), CNZ4(3.0, 3.0)
    oldqns = AbelianNumbers('G', [qn₁, qn₂, qn₃, qn₄, qn₂, qn₃, qn₁, qn₄], [2, 1, 3, 4, 1, 1, 2, 3], :counts)
    newqns, permutation = sort(oldqns)
    @test newqns == AbelianNumbers('C', [qn₁, qn₄, qn₂, qn₃], [4, 7, 2, 4], :counts)
    @test permutation == [1, 2, 13, 14, 7, 8, 9, 10, 15, 16, 17, 3, 11, 4, 5, 6, 12]
end

@testset "findall" begin
    qn₁, qn₂, qn₃ = CNZ4(2.0, 3.0), CNZ4(1.0, 2.0), CNZ4(0.0, 0.0)
    qns = AbelianNumbers('C', [qn₂, qn₁], [2, 3], :counts)
    @test findall(qn₁, qns, :compression) == [2]
    @test findall(qn₂, qns, :compression) == [1]
    @test findall(qn₃, qns, :compression) == []
    @test findall(qn₁, qns, :expansion) == [3, 4, 5]
    @test findall(qn₂, qns, :expansion) == [1, 2]
    @test findall(qn₃, qns, :expansion) == []

    qns = AbelianNumbers('G', [qn₁, qn₂, qn₁], [2, 2, 1], :counts)
    @test findall(qn₁, qns, :compression) == [1, 3]
    @test findall(qn₂, qns, :compression) == [2]
    @test findall(qn₃, qns, :compression) == []
    @test findall(qn₁, qns, :expansion) == [1, 2, 5]
    @test findall(qn₂, qns, :expansion) == [3, 4]
    @test findall(qn₃, qns, :expansion) == []
end

@testset "filter" begin
    qn₁, qn₂ = CNZ4(1.0, 2.0), CNZ4(3.0, 0.0)
    qns = AbelianNumbers('G', [qn₁, qn₂, qn₁, qn₂], [1, 2, 3, 4], :counts)
    @test filter((qn₁, qn₂), qns) == AbelianNumbers('G', [qn₁, qn₂, qn₁, qn₂], [1, 2, 3, 4], :counts)
    @test filter(qn₂, qns) == AbelianNumbers('G', [qn₂, qn₂], [2, 4], :counts)
end

@testset "prod" begin
    qn₁, qn₂ = CNZ4(1.0, 1.0), CNZ4(2.0, 3.0)
    qns₁ = AbelianNumbers('U', [qn₂, qn₁], [4, 5], :counts)
    qns₂ = AbelianNumbers('U', [qn₁, qn₂], [2, 3], :counts)
    qns, records = prod(qns₁, qns₂, signs=(+1, -1))
    @test qns == AbelianNumbers('C', [qn₁-qn₂, zero(CNZ4), qn₂-qn₁], [15, 22, 8], :counts)
    @test records == Dict(  (qn₁-qn₂) =>  Dict((qn₁, qn₂)=>1:15),
                            zero(CNZ4) => Dict((qn₂, qn₂)=>1:12, (qn₁, qn₁)=>13:22),
                            (qn₂-qn₁) =>  Dict((qn₂, qn₁)=>1:8)
                            )
end

@testset "expand" begin
    qn = CNZ4(3.0, 2.0)
    qns = AbelianNumbers(qn, 4)
    @test expand(qns, :indices) == [1, 1, 1, 1]
    @test expand(qns, :contents) == [qn, qn, qn, qn]
end

@testset "decompose" begin
    qn₁, qn₂ = CNZ4(0.0, 0.0), CNZ4(1.0, 1.0)
    qns = AbelianNumbers('U', [qn₁, qn₂], [0, 1, 2], :indptr)
    target = CNZ4(2.0, 2.0)
    result = Set(((1, 1, 2, 2), (1, 2, 1, 2), (1, 2, 2, 1), (2, 1, 1, 2), (2, 1, 2, 1), (2, 2, 1, 1)))
    @test ⊆(Set(decompose(target, qns, qns, qns, qns; signs=(1, 1, 1, 1), nmax=10, method=:bruteforce)), result)
    @test ⊆(Set(decompose(target, qns, qns, qns, qns; signs=(1, 1, 1, 1), nmax=10, method=:montecarlo)), result)
end

@testset "permute" begin
    qn₁, qn₂, qn₃ = CNZ4(1.0, 2.0), CNZ4(3.0, 0.0), CNZ4(4.0, 1.0)
    qns = AbelianNumbers('G', [qn₁, qn₂, qn₃], [2, 3, 4], :counts)
    @test permute(qns, [3, 2, 1], :compression) == AbelianNumbers('G', [qn₃, qn₂, qn₁], [4, 3, 2], :counts)
    @test permute(qns, [4, 6, 9, 8], :expansion) == AbelianNumbers('G', [qn₂, qn₃, qn₃, qn₃], [1, 1, 1, 1], :counts)
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

    @test Momentum₁{10} |> periods == (10,)
    @test Momentum₂{10, 15} |> periods == (10, 15)
    @test Momentum₃{10, 15, 20} |> periods == (10, 15, 20)

    @test Momentum₁{10}(1) == Momentum₁{10}(11) == Momentum₁{10}(-9)
    @test Momentum₂{10}(1, 1) == Momentum₂{10}(11, 11) == Momentum₂{10}(-9, -9)
    @test Momentum₃{10}(1, 1, 1) == Momentum₃{10}(11, 11, 11) == Momentum₃{10}(-9, -9, -9)
    @test Momentum₂{10, 20}(1, 1) == Momentum₂{10, 20}(11, 21) == Momentum₂{10, 20}(-9, -19)
    @test Momentum₃{10, 20, 30}(1, 1, 1) == Momentum₃{10, 20, 30}(11, 21, 31) == Momentum₃{10, 20, 30}(-9, -19, -29)
end