using Test
using StaticArrays: SVector
using LinearAlgebra: dot
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Essentials.DegreesOfFreedom
using QuantumLattices.Essentials.Spatials: PID,Point,pidtype,rcoord,icoord
using QuantumLattices.Mathematics.AlgebraOverFields: ID
import QuantumLattices.Interfaces: dimension,decompose,update!,sequence,reset!
import QuantumLattices.Essentials.DegreesOfFreedom: script,optdefaultlatex,latexsuperscript,latexsubscript
import QuantumLattices.Mathematics.AlgebraOverFields: rawelement

struct DID <: IID nambu::Int end
Base.adjoint(sl::DID)=DID(3-sl.nambu)

struct DIndex{S} <: Index{PID{S},DID}
    scope::S
    site::Int
    nambu::Int
end
Base.fieldnames(::Type{<:DIndex})=(:scope,:site,:nambu)
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:DID}=DIndex{fieldtype(P,:scope)}
function Base.angle(id::OID{<:DIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    phase=  length(vectors)==1 ? 2pi*dot(decompose(id.icoord,vectors[1]),values) :
            length(vectors)==2 ? 2pi*dot(decompose(id.icoord,vectors[1],vectors[2]),values) :
            length(vectors)==3 ? 2pi*dot(decompose(id.icoord,vectors[1],vectors[2],vectors[3]),values) :
            error("angle error: not supported number of input basis vectors.")
    id.index.nambu==1 ? phase : -phase
end

script(oid::OID{<:DIndex},::Val{:site})=oid.index.site
script(oid::OID{<:DIndex},::Val{:nambu})=oid.index.nambu==2 ? "\\dagger" : ""

struct DFock <: Internal{DID}
    atom::Int
    nnambu::Int
end
dimension(f::DFock)=f.nnambu
Base.getindex(f::DFock,i::Int)=DID(i)

struct DOperator{V,I<:ID} <: Operator{V,I}
    value::V
    id::I
end
rawelement(::Type{<:DOperator})=DOperator
optdefaultlatex(::Type{<:DOperator})=LaTeX{(:nambu,),(:site,)}('d')

@testset "Index" begin
    index=DIndex(PID('S',4),DID(1))
    @test index|>pidtype==PID{Char}
    @test index|>typeof|>pidtype==PID{Char}
    @test index|>iidtype==DID
    @test index|>typeof|>iidtype==DID
    @test index|>pid==PID('S',4)
    @test index|>iid==DID(1)
    @test index|>adjoint==DIndex('S',4,2)
    @test union(PID,IID)==Index{PID,IID}
    @test union(index|>pidtype,index|>iidtype)==index|>typeof
end

@testset "IndexToTuple" begin
    index=DIndex(PID('S',4),DID(1))
    @test directindextotuple(index)==('S',4,1)
    filteredindextotuple=FilteredAttributes(DIndex)
    @test filteredindextotuple==FilteredAttributes(:scope,:site,:nambu)
    @test isequal(filteredindextotuple,FilteredAttributes(:scope,:site,:nambu))
    @test filteredindextotuple|>length==3
    @test filteredindextotuple|>typeof|>length==3
    @test filteredindextotuple(index)==('S',4,1)
    @test filter(attr->attr≠:scope,filteredindextotuple)(index)==(4,1)
    @test filter(attr->attr≠:nambu,filteredindextotuple)(index)==('S',4)
    @test filter(attr->attr∉(:site,:nambu),filteredindextotuple)(index)==('S',)
end

@testset "Internal" begin
    it=DFock(1,2)
    @test it|>eltype==DID
    @test it|>typeof|>eltype==DID
    @test it==deepcopy(it)
    @test isequal(it,deepcopy(it))
    @test it|>string=="DFock(atom=1,nnambu=2)"
    @test it|>collect==[DID(1),DID(2)]
end

@testset "IDFConfig" begin
    config=IDFConfig{DFock}(pid->DFock((pid.site-1)%2+1,2),[PID(1,1),PID(1,2)])
    @test convert(Dict,config)==Dict(PID(1,1)=>DFock(1,2),PID(1,2)=>DFock(2,2))
    reset!(config,(PID(2,1),PID(2,2)))
    @test convert(Dict,config)==Dict(PID(2,1)=>DFock(1,2),PID(2,2)=>DFock(2,2))
end

@testset "Table" begin
    by=filter(attr->attr≠:nambu,FilteredAttributes(DIndex))

    table=Table([DIndex(1,1,1),DIndex(1,1,2)])
    @test table==Dict(DIndex(1,1,1)=>1,DIndex(1,1,2)=>2)
    table=Table([DIndex(1,1,1),DIndex(1,1,2)],by)
    @test table==Dict(DIndex(1,1,1)=>1,DIndex(1,1,2)=>1)

    config=IDFConfig{DFock}(pid->DFock((pid.site-1)%2+1,2),[PID(1,1),PID(1,2)])
    inds1=(DIndex(PID(1,1),iid) for iid in DFock(1,2))|>collect
    inds2=(DIndex(PID(1,2),iid) for iid in DFock(2,2))|>collect
    @test Table(config)==Table([inds1;inds2])
    @test Table(config,by)==Table([inds1;inds2],by)
    @test Table(config)==union(Table(inds1),Table(inds2))
    @test Table(config,by)|>reverse==Dict(1=>Set([DIndex(1,1,1),DIndex(1,1,2)]),2=>Set([DIndex(1,2,1),DIndex(1,2,2)]))

    table=Table(config)
    @test reset!(empty(table),[inds1;inds2])==table
    @test reset!(empty(table),config)==table
end

@testset "OID" begin
    oid=OID(DIndex(1,1,1),rcoord=SVector(0.0,-0.0),icoord=SVector(0.0,0.0),seq=1)
    @test oid'==OID(DIndex(1,1,2),rcoord=SVector(0.0,0.0),icoord=SVector(0.0,0.0),seq=1)
    @test hash(oid,UInt(1))==hash(OID(DIndex(1,1,1),rcoord=SVector(0.0,0.0),icoord=SVector(0.0,0.0),seq=1),UInt(1))
    @test propertynames(ID{<:NTuple{2,OID}},true)==(:contents,:indexes,:rcoords,:icoords,:seqs)
    @test propertynames(ID{<:NTuple{2,OID}},false)==(:indexes,:rcoords,:icoords,:seqs)
    @test fieldnames(OID)==(:index,:rcoord,:icoord,:seq)
    @test string(oid)=="OID(DIndex(1,1,1),[0.0,0.0],[0.0,0.0],1)"
    @test ID(oid',oid)'==ID(oid',oid)
    @test isHermitian(ID(oid',oid))==true
    @test isHermitian(ID(oid,oid))==false
    @test oidtype(DID,Point{2,PID{Char}},Nothing,Val(true))==OID{DIndex{Char},SVector{2,Float},SVector{2,Float},Nothing}
    @test oidtype(DID,Point{2,PID{Char}},Table{DIndex{Char}},Val(true))==OID{DIndex{Char},SVector{2,Float},SVector{2,Float},Int}
    @test oidtype(DID,Point{2,PID{Char}},Nothing,Val(false))==OID{DIndex{Char},Nothing,Nothing,Nothing}
    @test oidtype(DID,Point{2,PID{Char}},Table{DIndex{Char}},Val(false))==OID{DIndex{Char},Nothing,Nothing,Int}
end

@testset "Operator" begin
    @test rawelement(Operator{V} where V)==Operator
    @test rawelement(DOperator{V} where V)==DOperator

    opt=DOperator(1.0im,(DIndex(1,2,2),DIndex(1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),icoords=(SVector(2.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    @test opt'==DOperator(-1.0im,(DIndex(1,1,2),DIndex(1,2,1)),rcoords=(SVector(0.0,0.0),SVector(1.0,0.0)),icoords=(SVector(0.0,0.0),SVector(2.0,0.0)),seqs=(1,2))
    @test isHermitian(opt)==false
    @test string(opt)=="DOperator(value=1.0im,id=ID(OID(DIndex(1,2,2),[1.0,0.0],[2.0,0.0],2),OID(DIndex(1,1,1),[0.0,0.0],[0.0,0.0],1)))"
    @test twist(opt,[[1.0,0.0],[0.0,1.0]],[0.1,0.0])≈replace(opt,value=1.0im*conj(exp(2im*pi*0.2)))
    @test sequence(opt)==(2,1)
    @test sequence(opt,Dict(DIndex(1,2,2)=>3,DIndex(1,1,1)=>4))==(3,4)

    opt=DOperator(1.0,(DIndex(1,1,2),DIndex(1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.5,0.5)),icoords=(SVector(1.0,1.0),SVector(1.0,1.0)),seqs=(1,1))
    @test opt'==opt
    @test isHermitian(opt)==true

    opt=DOperator(1.0,(DIndex(1,1,2),DIndex(1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.0,0.5)),icoords=(SVector(1.0,1.0),SVector(0.0,1.0)),seqs=(1,1))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)

    opt=DOperator(1.0,ID(OID(DIndex(1,1,2),SVector(0.5,0.0),SVector(1.0,0.0),1)))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)
end

@testset "Operators" begin
    opt1=DOperator(1.0im,(DIndex(1,2,2),DIndex(1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),icoords=(SVector(2.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    opt2=DOperator(1.0,(DIndex(1,1,2),DIndex(1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),icoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1))
    opts=Operators(opt1,opt2)
    @test opts'==Operators(opt1',opt2')
    @test opts'+opts==Operators(opt1,opt1',opt2*2)
    @test isHermitian(opts)==false
    @test isHermitian(opts'+opts)==true
end

@testset "LaTeX" begin
    latex=LaTeX{(:nambu,),(:site,)}('c')
    @test latexsuperscript(latex|>typeof)==(:nambu,)
    @test latexsubscript(latex|>typeof)==(:site,)

    oid=OID(DIndex('d',1,2),rcoord=SVector(0.0,0.0),icoord=SVector(1.0,0.0),seq=1)
    @test script(oid,Val(:rcoord))=="[0.0,0.0]"
    @test script(oid,Val(:icoord))=="[1.0,0.0]"
    @test script(oid,latex,Val(:B))=='c'
    @test script(oid,latex,Val(:SP))==("\\dagger",)
    @test script(oid,latex,Val(:SB))==(1,)
    @test repr(oid,latex)=="c^{\\dagger}_{1}"

    opt=DOperator(1.0,(DIndex('d',2,2),DIndex('d',1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),icoords=(SVector(2.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    @test repr(opt,LaTeX{(:nambu,),(:site,)}('c'))=="c^{\\dagger}_{2}c^{}_{1}"
    @test repr(opt,LaTeX{(:nambu,),(:site,)}())=="d^{\\dagger}_{2}d^{}_{1}"
    @test repr(opt,LaTeX{(:nambu,),(:rcoord,)}())=="d^{\\dagger}_{[1.0,0.0]}d^{}_{[0.0,0.0]}"
    @test repr(opt)=="d^{\\dagger}_{2}d^{}_{1}"

    opts=Operators(DOperator(1.0-1.0im,(DIndex('d',2,2),DIndex('d',1,1))),DOperator(-1.0,(DIndex('d',1,2),DIndex('d',1,1))))
    @test repr(opts,LaTeX{(:nambu,),(:site,)}('c'))=="(1.0-1.0im)c^{\\dagger}_{2}c^{}_{1}-c^{\\dagger}_{1}c^{}_{1}"

    opt=DOperator(:h,(DIndex('d',2,2),DIndex('d',1,1)))
    @test repr(opt)==":hd^{\\dagger}_{2}d^{}_{1}"
end

@testset "Boundary" begin
    opt=DOperator(4.5,(DIndex('a',1,2),DIndex('b',2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(1,2))
    bound=Boundary{(:θ₁,:θ₂)}([0.1,0.2],[[1.0,0.0],[0.0,1.0]])
    @test bound==deepcopy(bound) && isequal(bound,deepcopy(bound))
    @test angle(bound,opt)≈0.6pi
    @test bound(opt)≈replace(opt,value=4.5*exp(2im*pi*0.3))
    update!(bound,θ₁=0.3)
    @test bound(opt)≈replace(opt,value=4.5*exp(2im*pi*0.5))
    bound=Boundary()
    @test angle(bound,opt)==0
    @test bound(opt)==opt
    @test update!(bound)==bound
end
