using Hamiltonian.Essentials.DegreesOfFreedom
using Hamiltonian.Essentials.Spatials: PID
import Hamiltonian.Prerequisites.Interfaces: dimension
import Hamiltonian.Essentials.DegreesOfFreedom: defaultcenter,propercenters

struct SlID <: IID nambu::Int end
Base.adjoint(sl::SlID)=SlID(3-sl.nambu)

struct SlIndex{S} <: Index{PID{S},SlID}
    scope::S
    site::Int
    nambu::Int
end
Base.fieldnames(::Type{<:SlIndex})=(:scope,:site,:nambu)
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:SlID}=SlIndex{fieldtype(P,:scope)}

struct SlFock <: Internal{SlID}
    atom::Int
    nnambu::Int
end
dimension(f::SlFock)=f.nnambu
Base.getindex(f::SlFock,i::Int)=SlID(i)

@testset "Index" begin
    index=SlIndex(PID('S',4),SlID(1))
    @test index|>pidtype==PID{Char}
    @test index|>typeof|>pidtype==PID{Char}
    @test index|>iidtype==SlID
    @test index|>typeof|>iidtype==SlID
    @test index|>pid==PID('S',4)
    @test index|>iid==SlID(1)
    @test index|>adjoint==SlIndex('S',4,2)
    @test union(PID,IID)==Index{PID,IID}
    @test union(index|>pidtype,index|>iidtype)==index|>typeof
end

@testset "IndexToTuple" begin
    index=SlIndex(PID('S',4),SlID(1))
    @test directindextotuple(index)==('S',4,1)
    filteredindextotuple=FilteredAttributes(SlIndex)
    @test filteredindextotuple==FilteredAttributes(:scope,:site,:nambu)
    @test filteredindextotuple|>length==3
    @test filteredindextotuple|>typeof|>length==3
    @test filteredindextotuple(index)==('S',4,1)
    @test filter(attr->attr≠:scope,filteredindextotuple)(index)==(4,1)
    @test filter(attr->attr≠:nambu,filteredindextotuple)(index)==('S',4)
    @test filter(attr->attr∉(:site,:nambu),filteredindextotuple)(index)==('S',)
end

@testset "Internal" begin
    it=SlFock(1,2)
    @test it|>eltype==SlID
    @test it|>typeof|>eltype==SlID
    @test it==deepcopy(it)
    @test isequal(it,deepcopy(it))
    @test it|>string=="SlFock(atom=1,nnambu=2)"
    @test it|>collect==[SlID(1),SlID(2)]
end

@testset "IDFConfig" begin
    config=IDFConfig(pid->SlFock((pid.site-1)%2+1,2),SlFock,[PID(1,1),PID(1,2)])
    @test convert(Dict,config)==Dict(PID(1,1)=>SlFock(1,2),PID(1,2)=>SlFock(2,2))
    replace!(config,PID(2,1),PID(2,2))
    @test convert(Dict,config)==Dict(PID(2,1)=>SlFock(1,2),PID(2,2)=>SlFock(2,2))
end

@testset "Table" begin
    by=filter(attr->attr≠:nambu,FilteredAttributes(SlIndex))

    table=Table([SlIndex(1,1,1),SlIndex(1,1,2)])
    @test table==Dict(SlIndex(1,1,1)=>1,SlIndex(1,1,2)=>2)
    table=Table([SlIndex(1,1,1),SlIndex(1,1,2)],by=by)
    @test table==Dict(SlIndex(1,1,1)=>1,SlIndex(1,1,2)=>1)

    config=IDFConfig(pid->SlFock((pid.site-1)%2+1,2),SlFock,[PID(1,1),PID(1,2)])
    inds1=(SlIndex(PID(1,1),iid) for iid in SlFock(1,2))|>collect
    inds2=(SlIndex(PID(1,2),iid) for iid in SlFock(2,2))|>collect
    @test Table(config)==Table([inds1;inds2])
    @test Table(config,by=by)==Table([inds1;inds2],by=by)
    @test Table(config)==union(Table(inds1),Table(inds2))
    @test Table(config,by=by)|>reverse==Dict(1=>Set([SlIndex(1,1,1),SlIndex(1,1,2)]),2=>Set([SlIndex(1,2,1),SlIndex(1,2,2)]))
end

@testset "Subscript" begin
    sub=@subscript (x1,x2)=>(x1,4,4,x2) with x1<x2
    @test sub==deepcopy(sub)
    @test isequal(sub,deepcopy(sub))
    @test rank(sub)==rank(typeof(sub))==2
    @test dimension(sub)==dimension(typeof(sub))==4
    @test string(sub)=="(x1,x2)=>(x1,4,4,x2) with $(sub.identifier)"
    @test sub(Val('M'),1,2)==(1,4,4,2)
    @test sub(Val('C'),1,2)==true
    @test sub(Val('C'),2,1)==false

    sub=@subscript (x1,x2)=>(x1,x2,x1,x2)
    @test sub(Val('M'),1,2)==(1,2,1,2)
    @test sub(Val('C'),1,2)==true

    sub=Subscript{4}()
    @test rank(sub)==1 && dimension(sub)==4
    @test string(sub)=="*=>(*,*,*,*)"
    @test sub(Val('M'),2)==(2,2,2,2)
    @test sub(Val('C'),2)==true

    sub=Subscript((1,2,2,1))
    @test rank(sub)==0 && dimension(sub)==4
    @test string(sub)=="(1,2,2,1)"
    @test sub(Val('M'))==(1,2,2,1)
    @test sub(Val('C'))==true
end

@testset "Subscripts" begin
    sub1=@subscript (x1,)=>(x1,2) with x1<2
    sub2=@subscript (y1,y2)=>(y1,y2,y1,y2) with y1<y2
    subs=Subscripts(sub1,sub2)
    @test rank(subs)==rank(typeof(subs))==3
    @test rank(subs,1)==rank(typeof(subs),1)==1
    @test rank(subs,2)==rank(typeof(subs),2)==2
    @test dimension(subs)==dimension(typeof(subs))==6
    @test dimension(subs,1)==dimension(typeof(subs),1)==2
    @test dimension(subs,2)==dimension(typeof(subs),2)==4
    @test subs(Val('M'),(1,2,3))==(1,2,2,3,2,3)
    @test subs(Val('C'),(1,2,3))==true
    @test subs(Val('C'),(2,2,3))==false
    @test subs(Val('C'),(1,3,2))==false
    @test subs(Val('C'),(2,3,2))==false

    sub3=Subscript((6,6))
    @test subs==sub1*sub2
    @test subs*sub3==Subscripts(sub1,sub2,sub3)
    @test sub3*subs==Subscripts(sub3,sub1,sub2)
    @test subs*Subscripts(sub3)==Subscripts(sub1,sub2,sub3)
    @test expand(subs*sub3,(3,3,2,3,2,3,8,9))|>collect==[(1,2,1,2,1,2,6,6),(1,2,1,3,1,3,6,6),(1,2,2,3,2,3,6,6)]
end

@testset "centerrelated" begin
    @test defaultcenter(Coupling,2,2,Val(1))==1
    @test defaultcenter(Coupling,1,2,Val(1))==1
    @test propercenters(Coupling,('*','*'),Val(1))==(1,1)
    @test propercenters(Coupling,(1,2,3),Val(3))==(1,2,3)
end
