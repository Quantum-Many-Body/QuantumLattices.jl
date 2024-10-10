using Base.Iterators: product
using QuantumLattices: ⊕, ⊗, ⊠, decompose, dimension, rank
using QuantumLattices.QuantumNumbers

@testset "AbelianQuantumNumber" begin
    n = ℕ(1)
    @test values(n) == 1
    @test n==ℕ(1) && isequal(n, ℕ(1))
    @test n<ℕ(2) && isless(n, ℕ(2))
    @test periods(n) == periods(typeof(n)) == (Inf,)
    @test period(n) == period(typeof(n)) == Inf
    @test zero(n) == zero(typeof(n)) == ℕ(0)
    @test n⊗ℕ(2) == ℕ(3)
    @test +n==n && n+ℕ(2) == ℕ(3)
    @test -n==ℕ(-1) && n-ℕ(2)==ℕ(-1)
    @test inv(n)==inv(n, true)==-n && inv(n, false)==n
    @test hash(n, UInt(10)) == hash(1, UInt(10))
    @test string(n) == "ℕ(1)"
    @test iterate(n) == (n, nothing)
    @test isnothing(iterate(n, nothing))

    sz = 𝕊ᶻ(1/2)
    sp = n ⊠ sz
    @test values(sp) == (1, 1/2)
    @test sp == Abelian[ℕ ⊠ 𝕊ᶻ](1, 1/2) == Abelian[ℕ ⊠ 𝕊ᶻ]((1, 1/2)) == CompositeAbelianQuantumNumber(n, sz)
    @test hash(sp, UInt(1)) == hash((n.charge, sz.charge), UInt(1))
    @test string(sp) == "Abelian[ℕ ⊠ 𝕊ᶻ](1, 1/2)"
    @test zero(sp) == zero(typeof(sp)) == Abelian[ℕ ⊠ 𝕊ᶻ](0, 0)
    @test length(sp) == rank(sp) == rank(typeof(sp)) == 2
    @test sp[1]==sp[begin]==n && sp[2]==sp[end]==sz
    @test periods(sp) == periods(typeof(sp)) == (Inf, Inf)
    @test period(sp, 1) == period(sp, 2) == Inf
    @test +sp==sp && sp+Abelian[ℕ ⊠ 𝕊ᶻ](1, 1/2)==Abelian[ℕ ⊠ 𝕊ᶻ](2, 1)
    @test -sp==Abelian[ℕ ⊠ 𝕊ᶻ](-1, -1/2) && Abelian[ℕ ⊠ 𝕊ᶻ](2, 1)-Abelian[ℕ ⊠ 𝕊ᶻ](1, 1/2)==sp

    @test (ℕ(1) ⊠ 𝕊ᶻ(1/2)) ⊠ (ℕ(2) ⊠ 𝕊ᶻ(3/2)) == (ℕ(1) ⊠ 𝕊ᶻ(1/2) ⊠ ℕ(2)) ⊠ 𝕊ᶻ(3/2) == ℕ(1) ⊠ (𝕊ᶻ(1/2) ⊠ ℕ(2) ⊠ 𝕊ᶻ(3/2))
    @test (ℕ ⊠ 𝕊ᶻ) ⊠ (ℕ ⊠ 𝕊ᶻ) == (ℕ ⊠ 𝕊ᶻ ⊠ ℕ) ⊠ 𝕊ᶻ == ℕ ⊠ (𝕊ᶻ ⊠ ℕ ⊠ 𝕊ᶻ)

    @test ℤ₄(-4)==ℤ₄(0)==ℤ₄(4) && ℤ₄(-3)==ℤ₄(1)==ℤ₄(5) && ℤ₄(-2)==ℤ₄(2)==ℤ₄(6) && ℤ₄(-1)==ℤ₄(3)==ℤ₄(7)
    @test string(ℤ₂(0)) == "ℤ₂(0)" && string(ℤ₃(0)) == "ℤ₃(0)" && string(ℤ₄(0)) == "ℤ₄(0)" && string(ℤ{5}(0)) == "ℤ{5}(0)"

    @test Momentum₁{10} |> periods == (10,)
    @test Momentum₂{10, 15} |> periods == (10, 15)
    @test Momentum₃{10, 15, 20} |> periods == (10, 15, 20)

    @test Momentum₁{10}(1) == Momentum₁{10}(11) == Momentum₁{10}(-9)
    @test Momentum₂{10}(1, 1) == Momentum₂{10}(11, 11) == Momentum₂{10}(-9, -9)
    @test Momentum₃{10}(1, 1, 1) == Momentum₃{10}(11, 11, 11) == Momentum₃{10}(-9, -9, -9)
    @test Momentum₂{10, 20}(1, 1) == Momentum₂{10, 20}(11, 21) == Momentum₂{10, 20}(-9, -19)
    @test Momentum₃{10, 20, 30}(1, 1, 1) == Momentum₃{10, 20, 30}(11, 21, 31) == Momentum₃{10, 20, 30}(-9, -19, -29)

    @test string(Momentum₁{10}) == "Momentum₁{10}"
    @test string(Momentum₂{10, 15}) == "Momentum₂{10, 15}"
    @test string(Momentum₃{10, 15, 20}) == "Momentum₃{10, 15, 20}"
    @test string(Momentum₁{10}(2)) == "Momentum₁{10}(2)"
    @test string(Momentum₂{10, 20}(2, 3)) == "Momentum₂{10, 20}(2, 3)"
    @test string(Momentum₃{10, 20, 30}(2, 3, 4)) == "Momentum₃{10, 20, 30}(2, 3, 4)"

    @test Int(Momentum₁{10}(2)) == 3
    @test Int(Momentum₂{10, 20}(2, 3)) == 44
    @test Int(Momentum₃{10, 20, 30}(2, 3, 4)) == 1295
end

@testset "Momenta" begin
    momenta = Momenta(Momentum₂{2, 3})
    @test momenta == Momenta(Momentum₂{2, 3})
    @test momenta ≠ Momenta(Momentum₂{2, 4})
    @test isequal(momenta, Momenta(Momentum₂{2, 3}))
    @test !isequal(momenta, Momenta(Momentum₂{2, 4}))
    @test collect(momenta) == [Momentum₂{2, 3}(0, 0), Momentum₂{2, 3}(0, 1), Momentum₂{2, 3}(0, 2), Momentum₂{2, 3}(1, 0), Momentum₂{2, 3}(1, 1), Momentum₂{2, 3}(1, 2)]
    for momentum in momenta
        @test momenta[convert(CartesianIndex, momentum, momenta)] == momentum
    end 
end

@testset "regularize" begin
    quantumnumbers, dimensions, perm = regularize([ℕ(4), ℕ(2), ℕ(3), ℕ(1), ℕ(2)], [2, 3, 1, 4, 9])
    @test quantumnumbers == [ℕ(1), ℕ(2), ℕ(3), ℕ(4)]
    @test dimensions == [4, 12, 1, 2]
    @test perm == [4, 2, 5, 3, 1]
end

@testset "AbelianGradedSpace" begin
    qns = AbelianGradedSpace([ℕ(1), ℕ(2), ℕ(3), ℕ(4)], [4, 12, 1, 2]; ordercheck=true, duplicatecheck=true, degeneracycheck=true)
    @test string(qns) == "Graded{ℕ}(1=>4, 2=>12, 3=>1, 4=>2)"
    @test qns==Graded{ℕ}(1=>4, 2=>12, 3=>1, 4=>2)==Graded{ℕ}((1=>4, 2=>12, 3=>1, 4=>2))==Graded(ℕ(1)=>4, ℕ(2)=>12, ℕ(3)=>1, ℕ(4)=>2)==Graded((ℕ(1)=>4, ℕ(2)=>12, ℕ(3)=>1, ℕ(4)=>2))
    @test length(qns) == 4
    @test qns[1]==ℕ(1) && qns[2]==ℕ(2) && qns[3]==ℕ(3) && qns[4]==ℕ(4)
    @test qns[2:-1:1] == qns[[ℕ(2), ℕ(1)]] == Graded{ℕ}(1=>4, 2=>12)
    @test ℕ(1)∈qns && ℕ(2)∈qns && ℕ(3)∈qns && ℕ(4)∈qns && ℕ(5)∉qns
    @test dimension(qns)==19 && dimension(qns, 1)==dimension(qns, ℕ(1))==4 && dimension(qns, 2)==dimension(qns, ℕ(2))==12 && dimension(qns, 3)==dimension(qns, ℕ(3))==1 && dimension(qns, 4)==dimension(qns, ℕ(4))==2
    @test range(qns, 1)==range(qns, ℕ(1))==1:4 && range(qns, 2)==range(qns, ℕ(2))==5:16 && range(qns, 3)==range(qns, ℕ(3))==17:17 && range(qns, 4)==range(qns, ℕ(4))==18:19
    @test cumsum(qns, 0)==0 && cumsum(qns, 1)==cumsum(qns, ℕ(1))==4 && cumsum(qns, 2)==cumsum(qns, ℕ(2))==16 && cumsum(qns, 3)==cumsum(qns, ℕ(3))==17 && cumsum(qns, 4)==cumsum(qns, ℕ(4))==19
    @test collect(pairs(qns, dimension))==[ℕ(1)=>4, ℕ(2)=>12, ℕ(3)=>1, ℕ(4)=>2] && collect(pairs(qns, range))==[ℕ(1)=>1:4, ℕ(2)=>5:16, ℕ(3)=>17:17, ℕ(4)=>18:19]
    @test [findindex(i, qns, guess) for (i, guess) in zip(1:dimension(qns), [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4])] == [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4]

    qns = qns'
    @test string(qns)=="Graded{ℕ}(1=>4, 2=>12, 3=>1, 4=>2)'"
    @test qns[1]==ℕ(-1) && qns[2]==ℕ(-2) && qns[3]==ℕ(-3) && qns[4]==ℕ(-4)
    @test qns[2:-1:1] == qns[[ℕ(-2), ℕ(-1)]] == Graded{ℕ}(1=>4, 2=>12; dual=true)
    @test ℕ(-1)∈qns && ℕ(-2)∈qns && ℕ(-3)∈qns && ℕ(-4)∈qns && ℕ(-5)∉qns
    @test dimension(qns)==19 && dimension(qns, 1)==dimension(qns, ℕ(-1))==4 && dimension(qns, 2)==dimension(qns, ℕ(-2))==12 && dimension(qns, 3)==dimension(qns, ℕ(-3))==1 && dimension(qns, 4)==dimension(qns, ℕ(-4))==2
    @test range(qns, 1)==range(qns, ℕ(-1))==1:4 && range(qns, 2)==range(qns, ℕ(-2))==5:16 && range(qns, 3)==range(qns, ℕ(-3))==17:17 && range(qns, 4)==range(qns, ℕ(-4))==18:19
    @test cumsum(qns, 0)==0 && cumsum(qns, 1)==cumsum(qns, ℕ(-1))==4 && cumsum(qns, 2)==cumsum(qns, ℕ(-2))==16 && cumsum(qns, 3)==cumsum(qns, ℕ(-3))==17 && cumsum(qns, 4)==cumsum(qns, ℕ(-4))==19
    @test collect(pairs(qns, dimension))==[ℕ(-1)=>4, ℕ(-2)=>12, ℕ(-3)=>1, ℕ(-4)=>2] && collect(pairs(qns, range))==[ℕ(-1)=>1:4, ℕ(-2)=>5:16, ℕ(-3)=>17:17, ℕ(-4)=>18:19]
    @test [findindex(i, qns, guess) for (i, guess) in zip(1:dimension(qns), [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4])] == [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4]
end

@testset "AbelianGradedSpaceSum" begin
    qns₁, qns₂, qns₃ = Graded{ℕ}(1=>2, 2=>4, 4=>1), Graded{ℕ}(1=>1, 2=>4, 3=>1), Graded{ℕ}(1=>1, 2=>4, 4=>1)
    qns = AbelianGradedSpaceSum(qns₁, qns₂, qns₃)
    @test string(qns) == "Graded{ℕ}(1=>2, 2=>4, 4=>1) ⊕ Graded{ℕ}(1=>1, 2=>4, 3=>1) ⊕ Graded{ℕ}(1=>1, 2=>4, 4=>1)"
    @test rank(qns) == rank(typeof(qns)) == 3
    @test qns == qns₁ ⊕ qns₂ ⊕ qns₃ == (qns₁ ⊕ qns₂) ⊕ qns₃ == qns₁ ⊕ (qns₂ ⊕ qns₃) == (qns₁ ⊕ qns₂) ⊕ AbelianGradedSpaceSum(qns₃)
    @test dimension(qns) == 19
    @test [dimension(qns, i) for i = 1:length(qns)] == [2, 4, 1, 1, 4, 1, 1, 4, 1]
    @test [range(qns, i) for i = 1:length(qns)] == [1:2, 3:6, 7:7, 8:8, 9:12, 13:13, 14:14, 15:18, 19:19]
    @test collect(pairs(qns, dimension)) == [ℕ(1)=>2, ℕ(2)=>4, ℕ(4)=>1, ℕ(1)=>1, ℕ(2)=>4, ℕ(3)=>1, ℕ(1)=>1, ℕ(2)=>4, ℕ(4)=>1]
    @test collect(pairs(qns, range)) == [ℕ(1)=>1:2, ℕ(2)=>3:6, ℕ(4)=>7:7, ℕ(1)=>8:8, ℕ(2)=>9:12, ℕ(3)=>13:13, ℕ(1)=>14:14, ℕ(2)=>15:18, ℕ(4)=>19:19]
    @test decompose(qns; expand=false) == (Graded{ℕ}(1=>4, 2=>12, 3=>1, 4=>2), [1, 4, 7, 2, 5, 8, 6, 3, 9])
    @test decompose(qns; expand=true) == (Graded{ℕ}(1=>4, 2=>12, 3=>1, 4=>2), [1, 2, 8, 14, 3, 4, 5, 6, 9, 10, 11, 12, 15, 16, 17, 18, 13, 7, 19])
end

@testset "AbelianGradedSpaceProd" begin
    qns₁, qns₂, qns₃ = Graded{𝕊ᶻ}(-1/2=>1, 1/2=>2), Graded{𝕊ᶻ}(-1/2=>2, 1/2=>1), Graded{𝕊ᶻ}(-1/2=>2, 1/2=>2)
    qns = AbelianGradedSpaceProd(qns₁, qns₂, qns₃)
    @test string(qns) == "Graded{𝕊ᶻ}(-1/2=>1, 1/2=>2) ⊗ Graded{𝕊ᶻ}(-1/2=>2, 1/2=>1) ⊗ Graded{𝕊ᶻ}(-1/2=>2, 1/2=>2)"
    @test rank(qns) == rank(typeof(qns)) == 3
    @test qns == qns₁ ⊗ qns₂ ⊗ qns₃ == (qns₁ ⊗ qns₂) ⊗ qns₃ == qns₁ ⊗ (qns₂ ⊗ qns₃) == (qns₁ ⊗ qns₂) ⊗ AbelianGradedSpaceProd(qns₃)
    @test dimension(qns) == 36
    @test [dimension(qns, i) for i = 1:length(qns)] == [4, 4, 2, 2, 8, 8, 4, 4]
    @test [dimension(qns, i) for i in reverse.(reshape(collect(product(qns₃, qns₂, qns₁)), :))] == [4, 4, 2, 2, 8, 8, 4, 4]
    @test [range(qns, i) for i = 1:length(qns)] == [
        [1, 2, 5, 6], [3, 4, 7, 8], [9, 10], [11, 12], [13, 14, 17, 18, 25, 26, 29, 30], [15, 16, 19, 20, 27, 28, 31, 32], [21, 22, 33, 34], [23, 24, 35, 36]
    ]
    @test [range(qns, i) for i in reverse.(reshape(collect(product(qns₃, qns₂, qns₁)), :))] == [
        [1, 2, 5, 6], [3, 4, 7, 8], [9, 10], [11, 12], [13, 14, 17, 18, 25, 26, 29, 30], [15, 16, 19, 20, 27, 28, 31, 32], [21, 22, 33, 34], [23, 24, 35, 36]
    ]
    @test collect(pairs(qns, dimension)) == [𝕊ᶻ(-3/2) => 4, 𝕊ᶻ(-1/2) => 4, 𝕊ᶻ(-1/2) => 2, 𝕊ᶻ(1/2) => 2, 𝕊ᶻ(-1/2) => 8, 𝕊ᶻ(1/2) => 8, 𝕊ᶻ(1/2) => 4, 𝕊ᶻ(3/2) => 4]
    @test collect(pairs(qns, range)) == [
        𝕊ᶻ(-3/2) => [1, 2, 5, 6],
        𝕊ᶻ(-1/2) => [3, 4, 7, 8],
        𝕊ᶻ(-1/2) => [9, 10],
        𝕊ᶻ(1/2) => [11, 12],
        𝕊ᶻ(-1/2) => [13, 14, 17, 18, 25, 26, 29, 30],
        𝕊ᶻ(1/2) => [15, 16, 19, 20, 27, 28, 31, 32],
        𝕊ᶻ(1/2) => [21, 22, 33, 34],
        𝕊ᶻ(3/2) => [23, 24, 35, 36]
    ]
    @test decompose(qns; expand=false) == (Graded{𝕊ᶻ}(-3/2=>4, -1/2=>14, 1/2=>14, 3/2=>4), [1, 2, 3, 5, 4, 6, 7, 8])
    @test decompose(qns; expand=true) == (Graded{𝕊ᶻ}(-3/2=>4, -1/2=>14, 1/2=>14, 3/2=>4), [1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 13, 14, 17, 18, 25, 26, 29, 30, 11, 12, 15, 16, 19, 20, 27, 28, 31, 32, 21, 22, 33, 34, 23, 24, 35, 36])
    @test merge(qns) == (
        Graded{𝕊ᶻ}(-3/2=>4, -1/2=>14, 1/2=>14, 3/2=>4), 
        Dict(
            𝕊ᶻ(-3/2) => [(𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2))],
            𝕊ᶻ(-1/2) => [(𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(-1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(-1/2)), (𝕊ᶻ(1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2))],
            𝕊ᶻ(1/2) => [(𝕊ᶻ(-1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(-1/2))],
            𝕊ᶻ(3/2) => [(𝕊ᶻ(1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(1/2))]
            )
    )
    @test split(𝕊ᶻ(-3/2), qns; nmax=20) ⊆ split(𝕊ᶻ(-3/2), qns; nmax=Inf) == Set([(𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2))])
    @test split(𝕊ᶻ(-1/2), qns; nmax=20) ⊆ split(𝕊ᶻ(-1/2), qns; nmax=Inf) == Set([(𝕊ᶻ(1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2)), (𝕊ᶻ(-1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(-1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(-1/2))])
    @test split(𝕊ᶻ(1/2), qns; nmax=20) ⊆ split(𝕊ᶻ(1/2), qns; nmax=Inf) == Set([(𝕊ᶻ(-1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(1/2), 𝕊ᶻ(-1/2), 𝕊ᶻ(1/2)), (𝕊ᶻ(1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(-1/2))])
    @test split(𝕊ᶻ(3/2), qns; nmax=20) ⊆ split(𝕊ᶻ(3/2), qns; nmax=Inf) == Set([(𝕊ᶻ(1/2), 𝕊ᶻ(1/2), 𝕊ᶻ(1/2))])
end
