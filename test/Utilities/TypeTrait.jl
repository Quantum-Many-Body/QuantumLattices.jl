using Hamiltonian.Utilities.TypeTrait

struct WithTrait{F1,F2,F3}
    f1::F1
    f2::F2
    f3::F3
end
Base.:(==)(fc1::WithTrait,fc2::WithTrait) = ==(efficientoperations,fc1,fc2)
Base.isequal(fc1::WithTrait,fc2::WithTrait) = isequal(efficientoperations,fc1,fc2)
Base.:<(fc1::WithTrait,fc2::WithTrait) = <(efficientoperations,fc1,fc2)
Base.:isless(fc1::WithTrait,fc2::WithTrait) = isless(efficientoperations,fc1,fc2)
Base.replace(fc::WithTrait;kwargs...) = replace(efficientoperations,fc;kwargs...)

@testset "efficientoperations" begin
    @test ==(efficientoperations,(),())
    @test isequal(efficientoperations,(),())
    @test ==(efficientoperations,(1,2),(1,2,3))==false
    @test isequal(efficientoperations,(1,2),(1,2,3))==false

    fc1,fc2=WithTrait(1.0,2,3),WithTrait(1,2,3)
    @test fc1==fc2
    @test isequal(fc1,fc2)
    fc1,fc2=WithTrait(1.0,2,3),WithTrait(1,2,4.0)
    @test fc1<fc2
    @test isless(fc1,fc2)

    @test replace(WithTrait(1,2,3),f1='c')==WithTrait('c',2,3)
end
