using Hamiltonian.Essentials.DegreeOfFreedom
using Hamiltonian.Essentials.Spatial: PID

struct IFID <: IID
    orbital::Int
    spin::Int
    nambu::Int
end

struct IFIndex{S} <: Index{PID{S},IFID}
    scope::S
    site::Int
    orbital::Int
    spin::Int
    nambu::Int
end
Base.fieldnames(::Type{<:IFIndex})=(:scope,:site,:orbital,:spin,:nambu)
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:IFID}=IFIndex{fieldtype(P,:scope)}

@testset "Index" begin
    index=IFIndex(PID('S',1),IFID(2,3,4))
    @test index|>pidtype==PID{Char}
    @test index|>typeof|>pidtype==PID{Char}
    @test index|>iidtype==IFID
    @test index|>typeof|>iidtype==IFID
    @test index|>pid==PID('S',1)
    @test index|>iid==IFID(2,3,4)
    @test union(PID,IID)==Index{PID,IID}
    @test union(index|>pidtype,index|>iidtype)==index|>typeof
    @test convert(Tuple,index)==('S',1,2,3,4)
    @test convert(Tuple,index,mask=(:scope,))==(1,2,3,4)
    @test convert(Tuple,index,mask=(:spin,))==('S',1,2,4)
    @test convert(Tuple,index,mask=(:site,:nambu))==('S',2,3)
end

struct IFInternal <: Internal{IFID}
    atom::Int
    norbital::Int
    nspin::Int
    nnambu::Int
end
Base.length(f::IFInternal)=prod((f.norbital,f.nspin,f.nnambu))
function Base.iterate(f::IFInternal,state=1)
    if state>length(f)
        nothing
    else
        count=state-1
        nambu=count%f.nnambu+1
        count=count÷f.nnambu
        spin=count%f.nspin+1
        orbital=count÷f.nspin+1
        IFID(orbital,spin,nambu),state+1
    end
end

@testset "Internal" begin
    it=IFInternal(1,2,2,2)
    @test it|>eltype==IFID
    @test it|>typeof|>eltype==IFID
    @test it==deepcopy(it)
    @test isequal(it,deepcopy(it))
    @test it|>string=="IFInternal(atom=1,norbital=2,nspin=2,nnambu=2)"
    @test it|>collect==[IFID(1,1,1),IFID(1,1,2),IFID(1,2,1),IFID(1,2,2),IFID(2,1,1),IFID(2,1,2),IFID(2,2,1),IFID(2,2,2)]
end

@testset "IDFConfig" begin
    config=IDFConfig{IFInternal}(pid->IFInternal((pid.site-1)%2+1,1,2,2),directindextotuple,[PID(1,1),PID(1,2)])
    @test convert(Dict,config)==Dict(PID(1,1)=>IFInternal(1,1,2,2),PID(1,2)=>IFInternal(2,1,2,2))
    replace!(config,PID(2,1),PID(2,2))
    @test convert(Dict,config)==Dict(PID(2,1)=>IFInternal(1,1,2,2),PID(2,2)=>IFInternal(2,1,2,2))
end

@testset "Table" begin
    table=Table([IFIndex(1,2,1,1,1),IFIndex(1,2,1,1,2),IFIndex(1,1,1,1,1),IFIndex(1,1,1,1,2)])
    @test table==Dict(IFIndex(1,1,1,1,1)=>1,IFIndex(1,1,1,1,2)=>2,IFIndex(1,2,1,1,1)=>3,IFIndex(1,2,1,1,2)=>4)
    table=Table([IFIndex(1,2,1,1,1),IFIndex(1,2,1,1,2),IFIndex(1,1,1,1,1),IFIndex(1,1,1,1,2)],mask=(:nambu,))
    @test table==Dict(IFIndex(1,1,1,1,1)=>1,IFIndex(1,1,1,1,2)=>1,IFIndex(1,2,1,1,1)=>2,IFIndex(1,2,1,1,2)=>2)
    config=IDFConfig{IFInternal}(pid->IFInternal((pid.site-1)%2+1,1,1,2),directindextotuple,[PID(1,1),PID(1,2)])
    inds1=(IFIndex(PID(1,1),iid) for iid in IFInternal(1,1,1,2))|>collect
    inds2=(IFIndex(PID(1,2),iid) for iid in IFInternal(2,1,1,2))|>collect
    @test Table(config)==Table([inds1;inds2])
    @test Table(config,mask=(:nambu,))==Table([inds1;inds2],mask=(:nambu,))
    @test Table(config)==union(Table(inds1),Table(inds2))
    @test Table(config,mask=(:nambu,))|>reverse==Dict(1=>Set([IFIndex(1,1,1,1,1),IFIndex(1,1,1,1,2)]),2=>Set([IFIndex(1,2,1,1,1),IFIndex(1,2,1,1,2)]))
end
