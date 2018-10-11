using Test
using Hamiltonian.Utilities.NamedVector

@namedvector "PIDFF" (:scope,:site) Float64 Vector
@namedvector "PIDFV" (:scope,:site) Float64
@namedvector "PIDVF" (:scope,:site) (<:Real) Vector
@namedvector "PIDVV" (:scope,:site)

@testset "FixedDataFixedVector" begin
    pid=PIDFF(1.0,0.0)
    @test string(pid)=="PIDFF(1.0,0.0)"
    @test pid|>zero==PIDFF(0.0,0.0)
    @test pid[1]==pid.scope==1.0
    @test pid[2]==pid.site==0.0
    @test (pid[1]=2.0;pid[1]==2.0)
    @test (pid.site=3.0;pid.site==3.0)
    @test eltype(pid|>typeof,1)==Float64
    @test eltype(pid|>typeof,2)==Vector{Float64}
end

@testset "FixedDataVariedVector" begin
    pid=PIDFV(1.0,0.0)
    @test string(pid)=="PIDFV(1.0,0.0)"
    @test pid|>zero==PIDFV(0.0,0.0)
    @test pid[1]==pid.scope==1.0
    @test pid[2]==pid.site==0.0
    @test (pid[1]=2.0;pid[1]==2.0)
    @test (pid.site=3.0;pid.site==3.0)
    @test eltype(pid|>typeof,1)==Float64
    @test eltype(pid|>typeof,2)==Vector{Float64}
end

@testset "VariedDataFixedVector" begin
    pid=PIDVF(1.0,0.0)
    @test string(pid)=="PIDVF(1.0,0.0)"
    @test pid|>zero==PIDVF(0.0,0.0)
    @test pid[1]==pid.scope==1.0
    @test pid[2]==pid.site==0.0
    @test (pid[1]=2.0;pid[1]==2.0)
    @test (pid.site=3.0;pid.site==3.0)
    @test eltype(pid|>typeof,1)==Float64
    @test eltype(pid|>typeof,2)==Vector{Float64}
end

@testset "VariedDataVariedVector" begin
    pid=PIDVV(1.0,0.0)
    @test string(pid)=="PIDVV(1.0,0.0)"
    @test pid|>zero==PIDVV(0.0,0.0)
    @test pid[1]==pid.scope==1.0
    @test pid[2]==pid.site==0.0
    @test (pid[1]=2.0;pid[1]==2.0)
    @test (pid.site=3.0;pid.site==3.0)
    @test eltype(pid|>typeof,1)==Float64
    @test eltype(pid|>typeof,2)==Vector{Float64}
end

@testset "AllEqual" begin
    @test PIDFF(1.0,2.0)==PIDVF(1,2)==PIDFV(1.0,2.0)==PIDVV(1,2)
end
