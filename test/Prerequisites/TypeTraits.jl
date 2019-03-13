using Test
using QuantumLattices.Prerequisites.TypeTraits

struct WithTrait{F1,F2,F3}
    f1::F1
    f2::F2
    f3::F3
end
Base.:(==)(fc1::WithTrait,fc2::WithTrait) = ==(efficientoperations,fc1,fc2)
Base.isequal(fc1::WithTrait,fc2::WithTrait) = isequal(efficientoperations,fc1,fc2)
Base.:<(fc1::WithTrait,fc2::WithTrait) = <(efficientoperations,fc1,fc2)
Base.isless(fc1::WithTrait,fc2::WithTrait) = isless(efficientoperations,fc1,fc2)
Base.isapprox(fc1::WithTrait,fc2::WithTrait;atol::Real=10^-5,rtol::Real=10^-5) = isapprox(efficientoperations,Val((:f1,:f2)),fc1,fc2;atol=atol,rtol=rtol)
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

    fc1,fc2=WithTrait(1.0+10^-6,2,3),WithTrait(1,2,3)
    @test fc1≈fc2
    fc1,fc2=WithTrait(1.0,2-10^-6,3),WithTrait(1,2,3)
    @test fc1≈fc2
    fc1,fc2=WithTrait(1.0,2,3-10^-6),WithTrait(1,2,3)
    @test fc1≉ fc2

    @test replace(WithTrait(1,2,3),f1='c')==WithTrait('c',2,3)
end

@testset "forder/corder" begin
    dims=(2,2,2)
    finds=[(1,1,1),(2,1,1),(1,2,1),(2,2,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]
    cinds=[(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    @test collect(indtosub(dims,i,forder) for i=1:prod(dims))==finds
    @test collect(indtosub(dims,i,corder) for i=1:prod(dims))==cinds
    @test collect(subtoind(dims,inds,forder) for inds in finds)==collect(1:prod(dims))
    @test collect(subtoind(dims,inds,corder) for inds in cinds)==collect(1:prod(dims))
end
