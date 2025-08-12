using QuantumLattices: OneOrMore, ZeroOrMore

@testset "OneOrMore" begin
    @test isa(1, OneOrMore{Integer}) && isa((1,), OneOrMore{Integer}) && isa((1, 2), OneOrMore{Integer})
    @test OneOrMore(1)==(1,) && OneOrMore((1,))==(1,) && OneOrMore((1, 2))==(1, 2)
end

@testset "ZeroOrMore" begin
    @test isa((), ZeroOrMore{Integer}) && isa(1, ZeroOrMore{Integer}) && isa((1,), ZeroOrMore{Integer}) && isa((1, 2), ZeroOrMore{Integer})
    @test ZeroOrMore(())==() && ZeroOrMore(1)==(1,) && ZeroOrMore((1,))==(1,) && ZeroOrMore((1, 2))==(1, 2)
end