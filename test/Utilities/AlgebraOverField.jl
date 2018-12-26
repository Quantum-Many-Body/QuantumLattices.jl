using Hamiltonian.Utilities.AlgebraOverField
using Hamiltonian.Utilities: Float

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
    @test hash(cid,UInt(1))==hash(cid.contents,UInt(1))

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

@testset "VectorSpace" begin
    id1,id2,id3=SMPID(1,1),SMPID(1,2),SMPID(1,3)
    svs=SimpleVectorSpace(id1,id2)
    @test dimension(svs)==2
    @test svs==id1⊕id2
    @test svs⊕id3==id1⊕id2⊕id3
    @test id3⊕svs==id3⊕id1⊕id2
    @test (id2⊕id3)⊕svs==id2⊕id3⊕id1⊕id2
    @test VectorSpace(id1,id2)==svs

    svs1=SimpleVectorSpace(SMPID(1,1),SMPID(1,2))
    svs2=SimpleVectorSpace(SMPID(2,1),SMPID(2,2))
    svs3=SimpleVectorSpace(SMPID(3,1),SMPID(3,2))
    cvs=CompositeVectorSpace(svs1,svs2)
    @test dimension(cvs)==4
    @test cvs==svs1⊗svs2
    @test cvs⊗svs3==svs1⊗svs2⊗svs3
    @test svs3⊗cvs==svs3⊗svs1⊗svs2
    @test (svs2⊗svs3)⊗cvs==svs2⊗svs3⊗svs1⊗svs2
    @test VectorSpace(svs1,svs2)==cvs
end

struct BasicOperator{V,I,N} <: Element{V,I,N}
    value::V
    id::I
    BasicOperator(value::Number,id::ID)=new{value|>typeof,id|>typeof,id|>typeof|>rank}(value,id)
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
