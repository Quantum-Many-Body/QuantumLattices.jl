using Test
using QuantumLattices.Prerequisites.CompositeStructures

struct FT{S,N,T} <: CompositeNTuple{N,T}
    info::S
    contents::NTuple{N,T}
end

@testset "CompositeNTuple" begin
    t=FT("Info",(1,2,3,4))
    @test t|>length==4
    @test t|>typeof|>length==4
    @test t|>eltype==Int
    @test t|>typeof|>eltype==Int
    @test t==deepcopy(t)
    @test isequal(t,deepcopy(t))
    @test hash(t)==hash(deepcopy(t))
    @test t[1]==1
    @test t[1:3]==FT("Info",(1,2,3))
    @test collect(t)==[1,2,3,4]
    @test collect(t|>Iterators.reverse)==[4,3,2,1]
    @test t|>keys==(1,2,3,4)|>keys
    @test t|>values==(1,2,3,4)|>values
    @test t|>pairs==(1,2,3,4)|>pairs
    @test t|>reverse==FT("Info",(4,3,2,1))
    @test convert(Tuple,t)==(1,2,3,4)
end

struct FV{S,T} <: CompositeVector{T}
    info::S
    contents::Vector{T}
end

@testset "CompositeVector" begin
    v=FV("Info",[1,3,2,4])
    @test size(v)==(4,)
    @test size(v,1)==4
    @test length(v)==4
    @test v==deepcopy(v)
    @test isequal(v,deepcopy(v))
    @test v[2]==3
    @test v[CartesianIndex(2,)]==3
    @test v[[1,3,4]]==FV("Info",[1,2,4])
    @test v[1:3]==FV("Info",[1,3,2])
    @test (v[1]=10;v[2:3]=11:12;v==FV("Info",[10,11,12,4]))
    @test push!(v,1)==FV("Info",[10,11,12,4,1])
    @test push!(v,2,3)==FV("Info",[10,11,12,4,1,2,3])
    @test pushfirst!(v,20,21)==FV("Info",[20,21,10,11,12,4,1,2,3])
    @test insert!(v,2,30)==FV("Info",[20,30,21,10,11,12,4,1,2,3])
    @test append!(v,[0,-1])==FV("Info",[20,30,21,10,11,12,4,1,2,3,0,-1])
    @test prepend!(v,[8,9])==FV("Info",[8,9,20,30,21,10,11,12,4,1,2,3,0,-1])
    @test splice!(v,2)==9 && v==FV("Info",[8,20,30,21,10,11,12,4,1,2,3,0,-1])
    @test splice!(v,1,9)==8 && v==FV("Info",[9,20,30,21,10,11,12,4,1,2,3,0,-1])
    @test splice!(v,1:3)==FV("Info",[9,20,30]) && v==FV("Info",[21,10,11,12,4,1,2,3,0,-1])
    @test splice!(v,1:3,[8,7,6])==FV("Info",[21,10,11]) && v==FV("Info",[8,7,6,12,4,1,2,3,0,-1])
    @test deleteat!(v,4)==FV("Info",[8,7,6,4,1,2,3,0,-1])
    @test deleteat!(v,[1,2])==FV("Info",[6,4,1,2,3,0,-1])
    @test deleteat!(v,5:6)==FV("Info",[6,4,1,2,-1])
    @test pop!(v)==-1 && v==FV("Info",[6,4,1,2])
    @test popfirst!(v)==6 && v==FV("Info",[4,1,2])
    @test empty!(v)==FV("Info",Int[]) && v==FV("Info",Int[])
    v=FV("Info",[1,3,2,2,4])
    @test empty(v)==FV("Info",Int[])
    @test v|>collect==[1,3,2,2,4]
    @test v|>keys==[1,3,2,2,4]|>keys
    @test v|>values==[1,3,2,2,4]|>values
    @test v|>pairs==[1,3,2,2,4]|>pairs
    @test convert(Vector,v)==[1,3,2,2,4]
    @test v|>reverse==FV("Info",[4,2,2,3,1])
    @test sort(v)==FV("Info",[1,2,2,3,4]) && v==FV("Info",[1,3,2,2,4])
    @test sort!(v)==FV("Info",[1,2,2,3,4]) && v==FV("Info",[1,2,2,3,4])
    @test filter(x->x<=3,v)==FV("Info",[1,2,2,3]) && v==FV("Info",[1,2,2,3,4])
    @test filter!(x->x<=3,v)==FV("Info",[1,2,2,3]) && v==FV("Info",[1,2,2,3])
    @test findfirst(isequal(2),v)==2
    @test findlast(isequal(2),v)==3
    @test findall(isequal(2),v)==[2,3]
end

struct FD{S,P,I} <: CompositeDict{P,I}
    info::S
    contents::Dict{P,I}
end

@testset "CompositeDict" begin
    d=FD("Info",Dict("a"=>1,"b"=>2))
    @test d==deepcopy(d)
    @test isequal(d,deepcopy(d))
    @test isempty(d)==false
    @test length(d)==2
    @test haskey(d,"a")==true && haskey(d,"d")==false
    @test Pair("a",1) ∈ d && Pair("d",4) ∉ d
    @test get(d,"a",2)==1
    @test get(d,"d",4)==4
    @test get(()->4,d,"d")==4
    @test get!(d,"d",4)==4 && d==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4))
    @test get!(()->5,d,"d")==4 && d==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4))
    @test getkey(d,"e","e")=="e"
    @test d["d"]==4
    @test push!(d,Pair("e",4))==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4,"e"=>4)) && d==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4,"e"=>4))
    @test (d["d"]=4;d["e"]=5;d==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4,"e"=>5)))
    @test pop!(d)==Pair("e",5) && d==FD("Info",Dict("a"=>1,"b"=>2,"d"=>4))
    @test pop!(d,"a")==1 && d==FD("Info",Dict("b"=>2,"d"=>4))
    @test pop!(d,"a",1)==1 && d==FD("Info",Dict("b"=>2,"d"=>4))
    @test delete!(d,"b")==FD("Info",Dict("d"=>4)) && d==FD("Info",Dict("d"=>4))
    @test empty!(d)==FD("Info",Dict{String,Int}()) && d==FD("Info",Dict{String,Int}())
    d=FD("Info",Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4))
    @test merge(FD("Info",Dict("a"=>1,"b"=>2)),FD("Info",Dict("c"=>3,"d"=>4)))==d
    @test merge(+,FD("Info",Dict("a"=>1,"b"=>2,"c"=>1)),FD("Info",Dict("c"=>2,"d"=>4)))==d
    @test empty(d)==FD("Info",Dict{String,Int}())
    @test d|>collect==Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4)|>collect
    @test d|>keys|>collect==Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4)|>keys|>collect
    @test d|>values|>collect==Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4)|>values|>collect
    @test d|>pairs|>collect==Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4)|>pairs|>collect
    @test convert(Dict,d)==Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4)
    @test filter(p->p.second<=3,d)==FD("Info",Dict("a"=>1,"b"=>2,"c"=>3)) && d==FD("Info",Dict("a"=>1,"b"=>2,"c"=>3,"d"=>4))
    @test filter!(p->p.second<=3,d)==FD("Info",Dict("a"=>1,"b"=>2,"c"=>3)) && d==FD("Info",Dict("a"=>1,"b"=>2,"c"=>3))
end
