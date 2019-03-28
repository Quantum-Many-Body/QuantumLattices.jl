using Test
using Printf: @printf
using QuantumLattices.Mathematics.AlgebraOverFields
using QuantumLattices.Interfaces: dimension,rank,add!,sub!,mul!,div!
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.Combinatorics: Combinations
using QuantumLattices.Mathematics.VectorSpaces: DirectVectorSpace
import QuantumLattices.Interfaces: ⊗,⋅

struct AOFID{O<:Real,S<:Real} <: SimpleID
    orbital::O
    spin::S
end

@testset "ID" begin
    cid=AOFID(2,1)*AOFID(1,Inf)
    @test cid==ID(AOFID(2,1),AOFID(1,Inf))
    @test cid|>propertynames==(:contents,)
    @test cid.contents==(AOFID(2,1),AOFID(1,Inf))
    @test cid|>string=="ID(AOFID(2,1),AOFID(1,Inf))"
    @test cid|>eltype==AOFID{Int,<:Real}
    @test cid|>rank==2
    @test cid|>typeof|>rank==2

    sid=AOFID(1,1)
    @test sid*cid==ID(AOFID(1,1),AOFID(2,1),AOFID(1,Inf))
    @test cid*sid==ID(AOFID(2,1),AOFID(1,Inf),AOFID(1,1))
    @test cid*cid==ID(AOFID(2,1),AOFID(1,Inf),AOFID(2,1),AOFID(1,Inf))

    @test AOFID(1,2)<AOFID(1,Inf)
    @test ID(AOFID(2,1))<ID(AOFID(1,2),AOFID(1,Inf))
    @test isless(AOFID(1,2),AOFID(1,Inf))
    @test isless(ID(AOFID(2,1)),ID(AOFID(1,2),AOFID(1,Inf)))
    @test AOFID(2,1)*AOFID(1,Inf)<AOFID(2,1)*AOFID(2,Inf)

    cid=AOFID(2,1)*AOFID(3,4)
    @test cid==ID(AOFID,(2,3),(1,4))
    @test propertynames(cid|>typeof,false)==(:orbitals,:spins)
    @test propertynames(cid|>typeof,true)==(:contents,:orbitals,:spins)
    @test cid.orbitals==(2,3)
    @test cid.spins==(1,4)
end

@testset "IdSpace" begin
    sids=DirectVectorSpace{'F'}((AOFID(i,1) for i=1:4)...)
    idspace=IdSpace(Combinations,sids,Val((0,2,4)))
    for i=1:dimension(idspace)
        @test searchsortedfirst(idspace,idspace[i])==i
        @test findfirst(idspace[i],idspace)==i
    end
end

struct AOFOperator{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I}
    value::V
    id::I
end
⊗(m1::AOFOperator,m2::AOFOperator)=m1*m2
⋅(m1::AOFOperator,m2::AOFOperator)=m1*m2
Base.show(io::IO,opt::AOFOperator)=@printf io "AOFOperator(value=%s,id=%s)" opt.value opt.id
Base.repr(opt::AOFOperator)="AOFOperator($(opt.value),$(opt.id))"

@testset "Elements" begin
    opt=AOFOperator(1.0,ID(AOFID(1,1)))
    @test opt|>deepcopy==opt
    @test isequal(opt|>deepcopy,opt)
    @test isapprox(opt,replace(opt,value=opt.value+10^-6);atol=10^-5)
    @test opt|>valtype==Float
    @test opt|>typeof|>valtype==Float
    @test opt|>idtype==ID{Tuple{AOFID{Int,Int}}}
    @test opt|>typeof|>idtype==ID{Tuple{AOFID{Int,Int}}}
    @test opt|>rank==1
    @test opt|>typeof|>rank==1
    @test +opt==opt
    @test opt*2==2*opt==AOFOperator(2.0,ID(AOFID(1,1)))
    @test -opt==AOFOperator(-1.0,ID(AOFID(1,1)))
    @test opt/2==AOFOperator(0.5,ID(AOFID(1,1)))
    @test opt^2==opt*opt
    @test opt+empty(Elements)==empty(Elements)+opt==opt-empty(Elements)
    @test empty(Elements)-opt==-opt

    opt1=AOFOperator(1.0,ID(AOFID(1,1)))
    opt2=AOFOperator(2.0,ID(AOFID(1,2)))
    opts=Elements(opt1.id=>opt1,opt2.id=>opt2)
    @test string(opts)=="Elements with 2 entries:\n  AOFOperator(value=2.0,id=ID(AOFID(1,2)))\n  AOFOperator(value=1.0,id=ID(AOFID(1,1)))\n"
    @test repr(opts)=="Elements with 2 entries:\n  AOFOperator(2.0,ID(AOFID(1,2)))\n  AOFOperator(1.0,ID(AOFID(1,1)))"
    @test opts|>zero==Elements{opt|>idtype,opt|>typeof}()
    @test opts|>typeof|>zero==Elements{opt|>idtype,opt|>typeof}()
    @test add!(deepcopy(opts))==opts
    @test add!(deepcopy(opts),empty(Elements))==opts
    @test sub!(deepcopy(opts))==opts
    @test sub!(deepcopy(opts),empty(Elements))==opts
    @test mul!(deepcopy(opts),2.0)==opts*2
    @test div!(deepcopy(opts),2.0)==opts/2
    @test +opts==opts+empty(Elements)==empty(Elements)+opts==opts==opts-empty(Elements)
    @test opt1+opt2==opts
    @test opt1+opts+opt2==Elements(opt1*2,opt2*2)
    @test opts+opts==Elements(opt1*2,opt2*2)
    @test -opts==Elements(-opt1,-opt2)==empty(Elements)-opts
    @test opts|>zero==opt1-opt1==opts-opts
    @test opts-opt1==Elements(opt2)
    @test opt1-opts==Elements(-opt2)
    @test 2*opts==opts*2==Elements(opt1*2,opt2*2)
    @test opts/2==Elements(opt1/2,opt2/2)
    @test opts^2==opts*opts

    opt3=AOFOperator(3.0,ID(AOFID(1,3)))
    @test opt1*opt2==opt1⊗opt2==opt1⋅opt2===AOFOperator(2.0,ID(AOFID(1,1),AOFID(1,2)))
    @test opts*opt3==opts⊗opt3==opts⋅opt3==Elements(AOFOperator(3.0,ID(AOFID(1,1),AOFID(1,3))),AOFOperator(6.0,ID(AOFID(1,2),AOFID(1,3))))
    @test opt3*opts==opt3⊗opts==opt3⋅opts==Elements(AOFOperator(3.0,ID(AOFID(1,3),AOFID(1,1))),AOFOperator(6.0,ID(AOFID(1,3),AOFID(1,2))))
    @test opts*Elements(opt3)==opts⊗Elements(opt3)==opts⋅Elements(opt3)==Elements(AOFOperator(3.0,ID(AOFID(1,1),AOFID(1,3))),AOFOperator(6.0,ID(AOFID(1,2),AOFID(1,3))))
end
