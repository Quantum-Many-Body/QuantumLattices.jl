using QuantumLattices: OneOrMore

@testset "OneOrMore" begin
    @test isa(1, OneOrMore{Integer}) && isa((1,), OneOrMore{Integer}) && isa((1, 2), OneOrMore{Integer})
    @test OneOrMore(1)==(1,) && OneOrMore((1,))==(1,) && OneOrMore((1, 2))==(1, 2)
end