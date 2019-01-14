using Hamiltonian.Mathematics.AlgebraOverFields
using Hamiltonian.Prerequisites: Float
using Hamiltonian.Prerequisites.Interfaces: dimension
using Hamiltonian.Mathematics.Combinatorics: Combinations
using Hamiltonian.Mathematics.VectorSpaces: DirectVectorSpace

struct SMPID{O<:Real,S<:Real} <: SimpleID
    orbital::O
    spin::S
end

@testset "ID" begin
    cid=SMPID(2,1)⊗SMPID(1,Inf)
    @test cid==ID(SMPID(2,1),SMPID(1,Inf))
    @test cid|>propertynames==(:contents,)
    @test cid.contents==(SMPID(2,1),SMPID(1,Inf))
    @test cid|>string=="ID(SMPID(2,1),SMPID(1,Inf))"
    @test cid|>eltype==SMPID{Int,<:Real}
    @test cid|>rank==2
    @test cid|>typeof|>rank==2

    sid=SMPID(1,1)
    @test sid⊗cid==ID(SMPID(1,1),SMPID(2,1),SMPID(1,Inf))
    @test cid⊗sid==ID(SMPID(2,1),SMPID(1,Inf),SMPID(1,1))
    @test cid⊗cid==ID(SMPID(2,1),SMPID(1,Inf),SMPID(2,1),SMPID(1,Inf))

    @test SMPID(1,2)<SMPID(1,Inf)
    @test ID(SMPID(2,1))<ID(SMPID(1,2),SMPID(1,Inf))
    @test isless(SMPID(1,2),SMPID(1,Inf))
    @test isless(ID(SMPID(2,1)),ID(SMPID(1,2),SMPID(1,Inf)))
    @test SMPID(2,1)⊗SMPID(1,Inf)<SMPID(2,1)⊗SMPID(2,Inf)

    cid=SMPID(2,1)⊗SMPID(3,4)
    @test cid==ID(SMPID,(2,3),(1,4))
    @test propertynames(cid|>typeof,false)==(:orbitals,:spins)
    @test propertynames(cid|>typeof,true)==(:contents,:orbitals,:spins)
    @test cid.orbitals==(2,3)
    @test cid.spins==(1,4)
end

@testset "IdSpace" begin
    sids=DirectVectorSpace{'F'}((SMPID(i,1) for i=1:4)...)
    idspace=IdSpace(Combinations,sids,Val((0,2,4)))
    for i=1:dimension(idspace)
        @test searchsortedfirst(idspace,idspace[i])==i
        @test findfirst(idspace[i],idspace)==i
    end
end

struct BasicOperator{V,I} <: Element{V,I}
    value::V
    id::I
end

@testset "Elements" begin
    opt=BasicOperator(1.0,ID(SMPID(1,1)))
    @test opt|>deepcopy==opt
    @test isequal(opt|>deepcopy,opt)
    @test opt|>valtype==Float
    @test opt|>typeof|>valtype==Float
    @test opt|>idtype==ID{1,SMPID{Int,Int}}
    @test opt|>typeof|>idtype==ID{1,SMPID{Int,Int}}
    @test opt|>rank==1
    @test opt|>typeof|>rank==1
    @test +opt==opt
    @test opt*2==2*opt==BasicOperator(2.0,ID(SMPID(1,1)))
    @test -opt==BasicOperator(-1.0,ID(SMPID(1,1)))
    @test opt/2==BasicOperator(0.5,ID(SMPID(1,1)))

    opt1=BasicOperator(1.0,ID(SMPID(1,1)))
    opt2=BasicOperator(2.0,ID(SMPID(1,2)))
    opts=Elements(opt1,opt2)
    @test opts|>zero==Elements{opt|>idtype,opt|>typeof}()
    @test opts|>typeof|>zero==Elements{opt|>idtype,opt|>typeof}()
    @test add!(deepcopy(opts))==opts
    @test sub!(deepcopy(opts))==opts
    @test +opts==opts
    @test opt1+opt2==opts
    @test opt1+opts+opt2==Elements(opt1*2,opt2*2)
    @test opts+opts==Elements(opt1*2,opt2*2)
    @test -opts==Elements(-opt1,-opt2)
    @test opts|>zero==opt1-opt1==opts-opts
    @test opts-opt1==Elements(opt2)
    @test opt1-opts==Elements(-opt2)
    @test 2*opts==opts*2==Elements(opt1*2,opt2*2)
    @test opts/2==Elements(opt1/2,opt2/2)

    opt3=BasicOperator(3.0,ID(SMPID(1,3)))
    @test opt1*opt2==BasicOperator(2.0,ID(SMPID(1,1),SMPID(1,2)))
    @test opts*opt3==Elements(BasicOperator(3.0,ID(SMPID(1,1),SMPID(1,3))),BasicOperator(6.0,ID(SMPID(1,2),SMPID(1,3))))
    @test opt3*opts==Elements(BasicOperator(3.0,ID(SMPID(1,3),SMPID(1,1))),BasicOperator(6.0,ID(SMPID(1,3),SMPID(1,2))))
    @test opts*Elements(opt3)==Elements(BasicOperator(3.0,ID(SMPID(1,1),SMPID(1,3))),BasicOperator(6.0,ID(SMPID(1,2),SMPID(1,3))))
end
