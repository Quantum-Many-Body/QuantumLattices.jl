using Hamiltonian.Utilities

@testset "indandsub" begin
    dims=(2,2,2)
    finds=[(1,1,1),(2,1,1),(1,2,1),(2,2,1),(1,1,2),(2,1,2),(1,2,2),(2,2,2)]
    cinds=[(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]
    @test collect(ind2sub(dims,i,forder) for i=1:prod(dims))==finds
    @test collect(ind2sub(dims,i,corder) for i=1:prod(dims))==cinds
    @test collect(sub2ind(dims,inds,forder) for inds in finds)==collect(1:prod(dims))
    @test collect(sub2ind(dims,inds,corder) for inds in cinds)==collect(1:prod(dims))
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
