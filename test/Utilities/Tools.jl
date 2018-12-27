using Hamiltonian.Utilities

@testset "indandsub" begin
    dims=(2,2,2)
    finds=[(1,1,1),(2,1,1),(1,2,1),(2,2,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]
    cinds=[(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    @test collect(indtosub(dims,i,forder) for i=1:prod(dims))==finds
    @test collect(indtosub(dims,i,corder) for i=1:prod(dims))==cinds
    @test collect(subtoind(dims,inds,forder) for inds in finds)==collect(1:prod(dims))
    @test collect(subtoind(dims,inds,corder) for inds in cinds)==collect(1:prod(dims))
end

@testset "decimaltostr" begin
    @test decimaltostr(1)=="1"
    @test decimaltostr(10^6)=="1000000"
    @test decimaltostr(1//7)=="1//7"
    @test decimaltostr(1.0)=="1.0"
    @test decimaltostr(10.0^6)=="1.0e+06"
    @test decimaltostr(10^-5)=="1.0e-05"
    @test decimaltostr(1/7)=="0.14286"
    @test decimaltostr(1/7,8)=="0.14285714"
    @test decimaltostr(0im)=="0.0"
    @test decimaltostr(1//7+0im)=="1//7"
    @test decimaltostr(1.0im)=="1.0im"
    @test decimaltostr(-1.0im)=="-1.0im"
    @test decimaltostr(0.1+0.12im)=="0.1+0.12im"
    @test decimaltostr(0.1-0.12im)=="0.1-0.12im"
end

@testset "ordinal" begin
    @test ordinal(1)=="1st"
    @test ordinal(2)=="2nd"
    @test ordinal(3)=="3rd"
    @test ordinal(4)=="4th"
    @test ordinal(5)=="5th"
end

struct ForTest{F1,F2,F3}
    f1::F1
    f2::F2
    f3::F3
end
Base.:(==)(fc1::ForTest,fc2::ForTest) = ==(efficientoperations,fc1,fc2)
Base.isequal(fc1::ForTest,fc2::ForTest) = isequal(efficientoperations,fc1,fc2)
Base.:<(fc1::ForTest,fc2::ForTest) = <(efficientoperations,fc1,fc2)
Base.:isless(fc1::ForTest,fc2::ForTest) = isless(efficientoperations,fc1,fc2)
Base.replace(fc::ForTest;kwargs...) = replace(efficientoperations,fc;kwargs...)

@testset "efficientoperations" begin
    @test ==(efficientoperations,(),())
    @test isequal(efficientoperations,(),())
    @test ==(efficientoperations,(1,2),(1,2,3))==false
    @test isequal(efficientoperations,(1,2),(1,2,3))==false

    fc1,fc2=ForTest(1.0,2,3),ForTest(1,2,3)
    @test fc1==fc2
    @test isequal(fc1,fc2)
    fc1,fc2=ForTest(1.0,2,3),ForTest(1,2,4.0)
    @test fc1<fc2
    @test isless(fc1,fc2)

    @test replace(ForTest(1,2,3),f1='c')==ForTest('c',2,3)
end

@testset "delta" begin
    @test delta(1,2)==0
    @test delta(1,1)==1
end
