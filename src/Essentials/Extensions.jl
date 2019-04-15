module Extensions

using ..Terms: Coupling,Couplings,subscriptexpr,Subscript
using ..FockPackage: CREATION,ANNIHILATION,FockCoupling
using ..SpinPackage: SpinCoupling
using ...Interfaces: ⊕,⊗,⋅

export @couplings,@fc_str,@sc_str
export @σ⁰_str,@σˣ_str,@σʸ_str,@σᶻ_str,@σ⁺_str,@σ⁻_str
export @heisenberg_str,@ising_str,@gamma_str,@sˣ_str,@sʸ_str,@sᶻ_str

"""
    @couplings cps -> Couplings

Convert an expression/literal to a set of couplings.
"""
macro couplings(cps)
    result=eval(cps)
    return isa(result,Coupling) ? Couplings(result) : isa(result,Couplings) ? result : error("@couplings error: inputs contain strangers that are not coupling/couplings.")
end

function cpcenters(str::AbstractString)
    @assert str[1]=='(' && str[end]==')' "cpcenters error: wrong input pattern."
    return Tuple(parse(Int,center) for center in split(str[2:end-1],','))
end

"""
    fc"..." -> FockCoupling

Construct a FockCoupling from a literal string.
"""
macro fc_str(str)
    ps=split(str," with ")
    conditions=length(ps)==2 ? fcconditions(ps[2]) : nothing
    ps=split(ps[1],' ')
    coeff=eval(Meta.parse(ps[1]))
    if ps[2][1]=='{' && ps[2][3]=='}'
        N=parse(Int,ps[2][2])
        return FockCoupling{N}(coeff)
    else
        attrpairs=[]
        components=split(ps[2],'@')
        centers=length(components)==2 ? cpcenters(components[2]) : nothing
        push!(attrpairs,:centers=>centers)
        N=centers===nothing ? nothing : length(centers)
        count=0
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=fccomponent(component)
                N===nothing && (N=length(attrvalue))
                if isa(attrvalue,Expr)
                    @assert attrname==:orbitals || attrname==:spins "@fc_str error: wrong input pattern."
                    @assert N==length(attrvalue.args) "@fc_str error: dismatched ranks."
                    count=count+1
                    condition=conditions===nothing ? :nothing : conditions[count]
                    push!(attrpairs,attrname=>eval(subscriptexpr(attrvalue,condition)))
                else
                    @assert N==length(attrvalue) "@fc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>attrvalue)
                end
            end
        end
        return FockCoupling{N}(coeff;attrpairs...)
    end
end
function fcconditions(str::AbstractString)
    conditions=Meta.parse(str)
    conditions=conditions.head==:tuple ? conditions.args : [conditions]
    return [condition=="*" ? :nothing : condition for condition in conditions]
end
function fccomponent(str::AbstractString)
    @assert str[3]=='(' && str[end]==')' "fccomponent error: wrong input pattern."
    attrname=str[1:2]=="sl" ? :atoms : str[1:2]=="ob" ? :orbitals : str[1:2]=="sp" ? :spins : str[1:2]=="ph" ? :nambus : error("fccomponent error: wrong input pattern.")
    expr=Meta.parse(str[3:end])
    attrvalue=isa(expr,Expr) ? (all(isa(arg,Int) for arg in expr.args) ? Tuple(expr.args) : expr) : (expr,)
    return attrname=>attrvalue
end

σᵅsplit(str::AbstractString)=(ps=split(str,'@'); length(ps)==1 ? (ps[1],nothing) : length(ps)==2 ? (ps[1],cpcenters(ps[2])) : "σᵅsplit error: wrong input pattern.")
σᵅname(mode::AbstractString)=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : mode=="ph" ? :nambus : error("σᵅname error: wrong input mode.")

"""
    σ⁰"sp"/σ⁰"sp@(c₁,c₂)" -> Couplings
    σ⁰"ob"/σ⁰"ob@(c₁,c₂)" -> Couplings
    σ⁰"sl"/σ⁰"sl@(c₁,c₂)" -> Couplings
    σ⁰"ph"/σ⁰"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁰_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σˣ"sp"/σˣ"sp@(c₁,c₂)" -> Couplings
    σˣ"ob"/σˣ"ob@(c₁,c₂)" -> Couplings
    σˣ"sl"/σˣ"sl@(c₁,c₂)" -> Couplings
    σˣ"ph"/σˣ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σˣ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σʸ"sp"/σʸ"sp@(c₁,c₂)" -> Couplings
    σʸ"ob"/σʸ"ob@(c₁,c₂)" -> Couplings
    σʸ"sl"/σʸ"sl@(c₁,c₂)" -> Couplings
    σʸ"ph"/σʸ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σʸ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(-1im;attrname=>attrval2,:centers=>centers)
end

"""
    σᶻ"sp"/σᶻ"sp@(c₁,c₂)" -> Couplings
    σᶻ"ob"/σᶻ"ob@(c₁,c₂)" -> Couplings
    σᶻ"sl"/σᶻ"sl@(c₁,c₂)" -> Couplings
    σᶻ"ph"/σᶻ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σᶻ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σ⁺"sp"/σ⁺"sp@(c₁,c₂)" -> Couplings
    σ⁺"ob"/σ⁺"ob@(c₁,c₂)" -> Couplings
    σ⁺"sl"/σ⁺"sl@(c₁,c₂)" -> Couplings
    σ⁺"ph"/σ⁺"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁺_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    σ⁻"sp"/σ⁻"sp@(c₁,c₂)" -> Couplings
    σ⁻"ob"/σ⁻"ob@(c₁,c₂)" -> Couplings
    σ⁻"sl"/σ⁻"sl@(c₁,c₂)" -> Couplings
    σ⁻"ph"/σ⁻"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁻_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    sc"..." -> SpinCoupling

Construct a SpinCoupling from a literal string.
"""
macro sc_str(str::String)
    ps=split(str," with ")
    condition=length(ps)==2 ? Meta.parse(ps[2]) : :nothing
    ps=split(ps[1],' ')
    coeff=eval(Meta.parse(ps[1]))
    tags=Tuple(replace(ps[2],"S"=>""))
    attrpairs=Any[:tags=>tags]
    if length(ps)==3
        components=split(ps[3],'@')
        centers=length(components)==2 ? cpcenters(components[2]) : nothing
        centers===nothing || @assert length(centers)==length(tags) "@sc_str error: dismatched ranks."
        push!(attrpairs,:centers=>centers)
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=sccomponent(component)
                if isa(attrvalue,Expr)
                    @assert attrname==:orbitals "@sc_str error: wrong input pattern."
                    @assert length(tags)==length(attrvalue.args) "@sc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>eval(subscriptexpr(attrvalue,condition)))
                else
                    @assert length(tags)==length(attrvalue) "@sc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>attrvalue)
                end
            end
        end
    end
    return SpinCoupling{length(tags)}(coeff;attrpairs...)
end
function sccomponent(str::AbstractString)
    @assert str[3]=='(' && str[end]==')' "sccomponent error: wrong input pattern."
    attrname=str[1:2]=="sl" ? :atoms : str[1:2]=="ob" ? :orbitals : error("sccomponent error: wrong input pattern.")
    expr=Meta.parse(str[3:end])
    attrvalue=isa(expr,Expr) ? (all(isa(arg,Int) for arg in expr.args) ? Tuple(expr.args) : expr) : (expr,)
    return attrname=>attrvalue
end

function scpairs(str::AbstractString,::Val{R}) where R
    attrpairs=Pair{Symbol,NTuple{R,Int}}[]
    if length(str)>0
        components=split(str,'@')
        length(components)==2 && push!(attrpairs,:centers=>cpcenters(components[2]))
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=sccomponent(component)
                @assert isa(attrvalue,NTuple{R,Int}) "scpairs error: wrong input pattern."
                push!(attrpairs,attrname=>attrvalue)
            end
        end
    end
    return attrpairs
end

"""
    heisenberg"sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Heisenberg couplings.
"""
macro heisenberg_str(str::String)
    attrpairs=scpairs(str,Val(2))
    sc1=SpinCoupling{2}(0.5;:tags=>('+','-'),attrpairs...)
    sc2=SpinCoupling{2}(0.5;:tags=>('-','+'),attrpairs...)
    sc3=SpinCoupling{2}(1.0;:tags=>('z','z'),attrpairs...)
    return Couplings(sc1,sc2,sc3)
end

"""
    ising"x sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    ising"y sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    ising"z sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Ising couplings.
"""
macro ising_str(str::String)
    @assert str[1] in ('x','y','z') "@ising_str error: wrong input pattern."
    attrpairs=length(str)>1 ? (@assert str[2]==' ' "@ising_str error: wrong input pattern."; scpairs(str[3:end],Val(2))) : Pair{Symbol,NTuple{2,Int}}[]
    return Couplings(SpinCoupling{2}(1.0;:tags=>(str[1],str[1]),attrpairs...))
end

"""
    gamma"x sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    gamma"y sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    gamma"z sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Gamma couplings.
"""
macro gamma_str(str::String)
    @assert str[1] in ('x','y','z') "@gamma_str error: wrong input pattern."
    t1,t2=str[1]=='x' ? ('y','z') : str[1]=='y' ? ('z','x') : ('x','y')
    attrpairs=length(str)>1 ? (@assert str[2]==' ' "@gamma_str error: wrong input pattern."; scpairs(str[3:end],Val(2))) : Pair{Symbol,NTuple{2,Int}}[]
    sc1=SpinCoupling{2}(1.0;:tags=>(t1,t2),attrpairs...)
    sc2=SpinCoupling{2}(1.0;:tags=>(t2,t1),attrpairs...)
    return Couplings(sc1,sc2)
end

"""
    sˣ"sl(a)⊗ob(o)" -> Couplings
    sʸ"sl(a)⊗ob(o)" -> Couplings
    sᶻ"sl(a)⊗ob(o)" -> Couplings

The single Sˣ/Sʸ/Sᶻ coupling.
"""
macro sˣ_str(str::String) Couplings(SpinCoupling{1}(1.0;:tags=>('x',),scpairs(str,Val(1))...)) end
macro sʸ_str(str::String) Couplings(SpinCoupling{1}(1.0;:tags=>('y',),scpairs(str,Val(1))...)) end
macro sᶻ_str(str::String) Couplings(SpinCoupling{1}(1.0;:tags=>('z',),scpairs(str,Val(1))...)) end

end
