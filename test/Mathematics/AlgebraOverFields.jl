using Test
using Printf: @printf
using QuantumLattices.Mathematics.AlgebraOverFields
using QuantumLattices.Interfaces: dimension,rank,add!,sub!,mul!,div!,sequence
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.Combinatorics: Combinations
using QuantumLattices.Mathematics.VectorSpaces: SimpleVectorSpace
import QuantumLattices.Interfaces: ⊗,⋅,permute

struct AOFID{O<:Real,S<:Real} <: SimpleID
    orbital::O
    nambu::S
end

@testset "ID" begin
    @test ID{SimpleID,1}|>rank==1
    @test ID|>rank==Any
    @test promote_type(ID{AOFID{Int,Int},1},ID{AOFID{Int,Int}})==ID{AOFID{Int,Int}}
    @test promote_type(ID{AOFID{Int,Int}},ID{AOFID{Int,Int},1})==ID{AOFID{Int,Int}}
    @test promote_type(ID{AOFID{Int,Int}},ID{AOFID{Int,Float}})==ID{AOFID{Int,<:Real}}
    @test promote_type(Tuple{},ID{AOFID{Int,Int}})==ID{AOFID{Int,Int}}
    @test promote_type(Tuple{},ID{AOFID{Int,Int},2})==ID{AOFID{Int,Int}}

    cid=AOFID(2,1)*AOFID(1,Inf)
    @test cid==ID(AOFID(2,1),AOFID(1,Inf))
    @test cid|>string=="ID(AOFID(2,1),AOFID(1,Inf))"
    @test cid|>eltype==AOFID{Int,<:Real}
    @test cid|>rank==2
    @test cid|>typeof|>rank==2

    sid=AOFID(1,1)
    @test sid*cid==ID(AOFID(1,1),AOFID(2,1),AOFID(1,Inf))
    @test cid*sid==ID(AOFID(2,1),AOFID(1,Inf),AOFID(1,1))
    @test cid*cid==ID(AOFID(2,1),AOFID(1,Inf),AOFID(2,1),AOFID(1,Inf))

    @test deepcopy(ID(AOFID(1,1)))==ID(AOFID(1,1))
    @test idisless(ID(AOFID(1,2)),ID(AOFID(1,Inf)))
    @test idisless(ID(AOFID(2,1)),ID(AOFID(1,2),AOFID(1,Inf)))

    cid=AOFID(2,1)*AOFID(3,4.0)
    @test cid==ID(AOFID,(2,3),(1,4.0))
    @test cid|>typeof|>idpropertynames==(:orbitals,:nambus)
    @test cid.orbitals==(2,3)
    @test cid.nambus==(1,4.0)
end

struct AOFOperator{V<:Number,I<:ID{SimpleID}} <: Element{V,I}
    value::V
    id::I
end
⊗(m1::AOFOperator,m2::AOFOperator)=m1*m2
⋅(m1::AOFOperator,m2::AOFOperator)=m1*m2
Base.show(io::IO,opt::AOFOperator)=@printf io "AOFOperator(value=%s,id=%s)" opt.value opt.id
Base.repr(opt::AOFOperator)="AOFOperator($(opt.value),$(opt.id))"

function permute(::Type{<:AOFOperator},id1::AOFID,id2::AOFID,::Any=nothing)
    @assert id1≠id2 "permute error: permuted ids should not be equal to each other."
    if id1.nambu==3-id2.nambu && id1.orbital==id2.orbital
        if id1.nambu==2
            return (AOFOperator(1,ID()),AOFOperator(1,ID(id2,id1)))
        else
            return (AOFOperator(-1,ID()),AOFOperator(1,ID(id2,id1)))
        end
    else
        return (AOFOperator(1,ID(id2,id1)),)
    end
end

@testset "Elements" begin
    @test valtype(Element)==Any
    @test valtype(Element{Int})==Int
    @test idtype(Element{<:Number})==ID{SimpleID}
    @test idtype(Element{<:Number,ID{AOFID}})==ID{AOFID}
    @test scalartype(AOFOperator{Number})==AOFOperator{Number,Tuple{}}
    @test promote_type(AOFOperator{Int},AOFOperator)==AOFOperator
    @test promote_type(AOFOperator,AOFOperator{Int})==AOFOperator
    @test promote_type(AOFOperator{Int},AOFOperator{Float})==AOFOperator{Float}
    @test promote_type(AOFOperator{Int,ID{AOFID{Int,Int},2}},AOFOperator{Float,ID{AOFID{Int,Int},2}})==AOFOperator{Float,ID{AOFID{Int,Int},2}}
    @test promote_type(AOFOperator{Int,ID{AOFID}},AOFOperator{Float,ID{AOFID}})==AOFOperator{Float,<:ID{AOFID}}

    opt=AOFOperator(1.0,ID(AOFID(1,1)))
    @test opt|>deepcopy==opt
    @test isequal(opt|>deepcopy,opt)
    @test isapprox(opt,replace(opt,value=opt.value+10^-6);atol=10^-5)
    @test opt|>valtype==opt|>typeof|>valtype==Float
    @test opt|>idtype==opt|>typeof|>idtype==ID{AOFID{Int,Int},1}
    @test opt|>rank==opt|>typeof|>rank==1
    @test opt|>typeof|>one==AOFOperator(1.0,ID())
    @test +opt==opt
    @test opt+zero(Elements)==zero(Elements)+opt==opt-zero(Elements)
    @test -opt==AOFOperator(-1.0,ID(AOFID(1,1)))
    @test zero(Elements)-opt==-opt
    @test opt*2==2*opt==opt*AOFOperator(2,ID())==AOFOperator(2,ID())*opt==AOFOperator(2.0,ID(AOFID(1,1)))
    @test opt*zero(Elements)==zero(Elements)*opt==zero(Elements)==nothing
    @test opt/2==opt/AOFOperator(2,ID())==AOFOperator(0.5,ID(AOFID(1,1)))
    @test opt^2==opt*opt

    @test AOFOperator(1,ID())+1==1+AOFOperator(1,ID())==AOFOperator(1,ID())+AOFOperator(1,ID())==AOFOperator(2,ID())
    @test opt+1==1+opt==opt+AOFOperator(1,ID())==AOFOperator(1,ID())+opt
    @test AOFOperator(1,ID())-2==1-AOFOperator(2,ID())==AOFOperator(1,ID())-AOFOperator(2,ID())==AOFOperator(-1,ID())
    @test AOFOperator(2,ID())*AOFOperator(3,ID())==AOFOperator(6,ID())

    opt1=AOFOperator(1.0,ID(AOFID(1,1)))
    opt2=AOFOperator(2.0,ID(AOFID(1,2)))
    opts=Elements(opt1.id=>opt1,opt2.id=>opt2)
    @test string(opts)=="Elements with 2 entries:\n  AOFOperator(value=2.0,id=ID(AOFID(1,2)))\n  AOFOperator(value=1.0,id=ID(AOFID(1,1)))\n"
    @test repr(opts)=="Elements with 2 entries:\n  AOFOperator(2.0,ID(AOFID(1,2)))\n  AOFOperator(1.0,ID(AOFID(1,1)))"
    @test opts|>deepcopy==opts
    @test opts|>zero==opts|>typeof|>zero==nothing
    @test add!(deepcopy(opts))==opts
    @test add!(deepcopy(opts),zero(Elements))==opts
    @test sub!(deepcopy(opts))==opts
    @test sub!(deepcopy(opts),zero(Elements))==opts
    @test mul!(deepcopy(opts),2.0)==mul!(deepcopy(opts),AOFOperator(2,ID()))==opts*2==opts*AOFOperator(2,ID())
    @test div!(deepcopy(opts),2.0)==div!(deepcopy(opts),AOFOperator(2,ID()))==opts/2==opts/AOFOperator(2,ID())
    @test +opts==opts+zero(Elements)==zero(Elements)+opts==opts==opts-zero(Elements)
    @test opt1+opt2==opts
    @test opt1+opts==opts+opt1==Elements(opt1*2,opt2)
    @test opts+opts==Elements(opt1*2,opt2*2)
    @test -opts==Elements(-opt1,-opt2)==zero(Elements)-opts
    @test 1.0-opts==AOFOperator(1.0,ID())-opts==-opts+1.0
    @test opts-1.0==opts-AOFOperator(1.0,ID())
    @test opt-1.0==opt-AOFOperator(1.0,ID())==opt+AOFOperator(-1.0,ID())
    @test 1.0-opt==AOFOperator(1.0,ID())-opt==1.0+(-opt)
    @test zero(Elements)==opt1-opt1==opts-opts
    @test opts-opts==zero(Elements)
    @test isequal(zero(Elements),opts-opts)
    @test isequal(opts-opts,zero(Elements))
    @test opts-opt1==Elements(opt2)
    @test opt1-opts==Elements(-opt2)
    @test 2*opts==opts*2==Elements(opt1*2,opt2*2)==AOFOperator(2,ID())*opts
    @test opts/2==Elements(opt1/2,opt2/2)
    @test opts^2==opts*opts
    @test opts*zero(Elements)==zero(Elements)*opts==nothing

    temp=Elements{ID{AOFID{Int,Int}},AOFOperator{Float}}(opts)
    @test add!(deepcopy(temp),1)==add!(deepcopy(temp),AOFOperator(1,ID()))==opts+1==1+opts==opts+AOFOperator(1,ID())==AOFOperator(1,ID())+opts
    @test sub!(deepcopy(temp),1)==sub!(deepcopy(temp),AOFOperator(1,ID()))==opts-1==opts-AOFOperator(1,ID())

    opt3=AOFOperator(3,ID(AOFID(1,3)))
    @test opt1*opt2==opt1⊗opt2==opt1⋅opt2==AOFOperator(2.0,ID(AOFID(1,1),AOFID(1,2)))
    @test opts*opt3==opts⊗opt3==opts⋅opt3==Elements(AOFOperator(3.0,ID(AOFID(1,1),AOFID(1,3))),AOFOperator(6.0,ID(AOFID(1,2),AOFID(1,3))))
    @test opt3*opts==opt3⊗opts==opt3⋅opts==Elements(AOFOperator(3.0,ID(AOFID(1,3),AOFID(1,1))),AOFOperator(6.0,ID(AOFID(1,3),AOFID(1,2))))
    @test opts*Elements(opt3)==opts⊗Elements(opt3)==opts⋅Elements(opt3)==Elements(AOFOperator(3.0,ID(AOFID(1,1),AOFID(1,3))),AOFOperator(6.0,ID(AOFID(1,2),AOFID(1,3))))

    table=Dict(AOFID(1,1)=>1,AOFID(1,2)=>2)
    opt=AOFOperator(2.0,ID(AOFID(1,1),AOFID(1,2)))
    @test sequence(opt,table)==(1,2)
    @test split(opt)==(2.0,AOFOperator(1.0,ID(AOFID(1,1))),AOFOperator(1.0,ID(AOFID(1,2))))
end

@testset "replace" begin
    id1,id2=AOFID(1,1),AOFID(1,2)
    opt=AOFOperator(1.5,ID(id1,id2))
    id3,id4=AOFID(2,1),AOFID(2,2)
    opt1=AOFOperator(2.0,ID(id3))+AOFOperator(3.0,ID(id4))
    opt2=AOFOperator(2.0,ID(id3))-AOFOperator(3.0,ID(id4))
    @test replace(opt,id1=>opt1,id2=>opt2)==1.5*opt1*opt2
    @test replace(Elements(opt),id2=>opt2,id1=>opt1)==1.5*opt1*opt2
end

@testset "permute" begin
    id1,id2=AOFID(1,1),AOFID(1,2)
    opt1=AOFOperator(1.5,ID(id1,id2))
    opt2=AOFOperator(1.5,ID(id2,id1))
    opt0=AOFOperator(1.5,ID())
    @test permute(opt1,Dict(id1=>1,id2=>2))==opt2-opt0
    @test permute(opt1+opt2,Dict(id1=>1,id2=>2))==2*opt2-opt0
    @test permute(opt2,Dict(id1=>2,id2=>1))==opt1+opt0
    @test permute(opt1+opt2,Dict(id1=>2,id2=>1))==2*opt1+opt0
end
