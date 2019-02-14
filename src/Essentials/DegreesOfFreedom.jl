module DegreesOfFreedom

using Printf: @printf
using ...Prerequisites.TypeTraits: efficientoperations
using ...Prerequisites.CompositeStructures: CompositeDict,CompositeTuple
using ...Mathematics.VectorSpaces: AbstractVectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements
using ...Prerequisites.TypeTraits: indtosub,corder
using ..Spatials: PID

import ..Spatials: pidtype
import ...Prerequisites.Interfaces: rank,dimension,expand

export rank,dimension,expand
export IID,Index,pidtype,pid,iidtype,iid
export IndexToTuple,DirectIndexToTuple,directindextotuple,FilteredAttributes
export Internal,IDFConfig,Table
export Subscript,Subscripts,@subscript
export Coupling,Couplings

"""
    IID

The id of an internal degree of freedom.
"""
abstract type IID <: SimpleID end

"""
    Internal

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: AbstractVectorSpace{I} end

"""
    show(io::IO,i::Internal)

Show an internal.
"""
Base.show(io::IO,i::Internal)=@printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i,name))" for name in i|>typeof|>fieldnames),",")

"""
    Index{P,I}

The complete index of a degree of freedom, which consist of the spatial part and the internal part.
"""
abstract type Index{P<:PID,I<:IID} <: SimpleID end

"""
    (INDEX::Type{<:Index})(pid::PID,iid::IID) -> INDEX

Get the corresponding index from a pid and an iid.
"""
(INDEX::Type{<:Index})(pid::PID,iid::IID)=INDEX(convert(Tuple,pid)...,convert(Tuple,iid)...)

"""
    pidtype(index::Index)
    pidtype(::Type{<:Index{P,I}}) where {P,I}

Get the type of the spatial part of an index.
"""
pidtype(index::Index)=index|>typeof|>pidtype
pidtype(::Type{<:Index{P,I}}) where {P,I}=P

"""
    pid(index::Index) -> PID

Get the spatial part of an index.
"""
@generated function pid(index::Index)
    exprs=[:(getfield(index,$i)) for i=1:fieldcount(index|>pidtype)]
    return :(pidtype(index).name.wrapper($(exprs...)))
end

"""
    iidtype(index::Index)
    iidtype(::Type{<:Index{P,I}}) where {P,I}

Get the type of the internal part of an index.
"""
iidtype(index::Index)=index|>typeof|>iidtype
iidtype(::Type{<:Index{P,I}}) where {P,I}=I

"""
    iid(index::Index) -> IID

Get the internal part of an index.
"""
@generated function iid(index::Index)
    exprs=[:(getfield(index,$i)) for i=fieldcount(index|>pidtype)+1:fieldcount(index)]
    return :(iidtype(index).name.wrapper($(exprs...)))
end

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
Base.adjoint(index::Index)=typeof(index)(index|>pid,index|>iid|>adjoint)

"""
    union(::Type{P},::Type{I}) where {P<:PID,I<:IID}

Combine a concrete `PID` type and a concrete `IID` type to a concrete `Index` type.
"""
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:IID}=Index{P,I}

"""
    IndexToTuple

The rules for converting an index to a tuple.

As a function, every instance should accept only one positional argument, i.e. the index to be converted to a tuple.
"""
abstract type IndexToTuple <:Function end

"""
    DirectIndexToTuple

Direct index to tuple.
"""
struct DirectIndexToTuple <: IndexToTuple end

"""
    (indextotuple::DirectIndexToTuple)(index::Index) -> Tuple

Convert an index to tuple directly.
"""
(indextotuple::DirectIndexToTuple)(index::Index)=convert(Tuple,index)

"""
    directindextotuple

Indicate that the conversion from an index to a tuple is direct.
"""
const directindextotuple=DirectIndexToTuple()

"""
    FilteredAttributes(::Type{I}) where I<:Index

A method that converts an arbitary index to a tuple, by iterating over the selected attributes in a specific order.
"""
struct FilteredAttributes{N} <: IndexToTuple
    attributes::NTuple{N,Symbol}
end
FilteredAttributes(attrs::Symbol...)=FilteredAttributes(attrs)
FilteredAttributes(::Type{I}) where I<:Index=FilteredAttributes(I|>fieldnames)

"""
    length(indextotuple::FilteredAttributes) -> Int
    length(::Type{<:FilteredAttributes{N}}) where N -> Int

Get the length of the filtered attributes.
"""
Base.length(indextotuple::FilteredAttributes)=indextotuple|>typeof|>length
Base.length(::Type{<:FilteredAttributes{N}}) where N=N

"""
    (indextotuple::FilteredAttributes)(index::Index) -> Tuple

Convert an index to tuple by the "filtered attributes" method.
"""
@generated function (indextotuple::FilteredAttributes)(index::Index)
    exprs=[:(getfield(index,indextotuple.attributes[$i])) for i=1:length(indextotuple)]
    return Expr(:tuple,exprs...)
end

"""
    filter(f::Function,indextotuple::FilteredAttributes)

Filter the attributes of a "filtered attributes" method.
"""
Base.filter(f::Function,indextotuple::FilteredAttributes)=FilteredAttributes(Tuple(attr for attr in indextotuple.attributes if f(attr)))

"""
    IDFConfig(map::Function,::Type{I},pids::AbstractVector{<:PID}=[]) where I<:Internal

Configuration of the internal degrees of freedom at a lattice.

`map` maps a `PID` to an `Internal`.
"""
struct IDFConfig{M<:Function,P<:PID,I<:Internal} <: CompositeDict{P,I}
    map::M
    contents::Dict{P,I}
end
function IDFConfig(map::Function,::Type{I},pids::AbstractVector{<:PID}=[]) where I<:Internal
    contents=Dict{pids|>eltype,I}()
    for pid in pids
        contents[pid]=map(pid)
    end
    IDFConfig(map,contents)
end

"""
    replace!(config::IDFConfig,pids::PID...) -> IDFConfig

Reset the idfconfig with new pids.
"""
function Base.replace!(config::IDFConfig,pids::PID...)
    empty!(config)
    for pid in pids
        config[pid]=config.map(pid)
    end
    config
end

"""
    Table{I<:Index} <: AbstractDict{I,Int}

Index-sequence table. Alias for `Dict{I<:Index,Int}`.
"""
const Table{I<:Index}=Dict{I,Int}

"""
    Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple) -> Table

Convert an sequence of indices to the corresponding index-sequence table.

The input indices will be converted to tuples by the `by` function with the duplicates removed. The resulting unique tuples are sorted, which determines the sequence of the input indices. Note that two indices have the same sequence if their converted tupels are equal to each other.
"""
function Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple)
    tuples=[by(index) for index in indices]
    permutation=sortperm(tuples,alg=Base.Sort.QuickSort)
    result=Table{indices|>eltype}()
    count=1
    for i=1:length(tuples)
        i>1 && tuples[permutation[i]]!=tuples[permutation[i-1]] && (count+=1)
        result[indices[permutation[i]]]=count
    end
    result
end

"""
    Table(config::IDFConfig;by::IndexToTuple=directindextotuple) -> Table

Get the index-sequence table of the whole internal Hilbert spaces at a lattice.
"""
function Table(config::IDFConfig;by::IndexToTuple=directindextotuple)
    result=union(config|>keytype,config|>valtype|>eltype)[]
    for (pid,internal) in config
        for iid in internal
            push!(result,(result|>eltype)(pid,iid))
        end
    end
    Table(result,by=by)
end

"""
    union(tables::Table...;by::IndexToTuple=directindextotuple) -> Table

Unite several index-sequence tables.

See [`Table`](@ref) for more details.
"""
function Base.union(tables::Table...;by::IndexToTuple=directindextotuple)
    indices=(tables|>eltype|>keytype)[]
    for table in tables
        for index in keys(table)
            push!(indices,index)
        end
    end
    Table(indices,by=by)
end

"""
    reverse(table::Table) -> Dict{Int,Set{<:Index}}

Convert an index-sequence table to a sequence-indices table.

Since different indices may correspond to the same sequence, the reverse is a one-to-many map.
"""
function Base.reverse(table::Table)
    result=Dict{Int,Set{table|>keytype}}()
    for (index,seq) in table
        haskey(result,seq) || (result[seq]=Set{table|>keytype}())
        push!(result[seq],index)
    end
    result
end

const wildcard='*'
const constant=':'
"""
    Subscript(  ipattern::NTuple{N1,Any},
                opattern::NTuple{N2,Any},
                mapping::Union{Function,Nothing}=nothing,
                constrain::Union{Function,Nothing}=nothing,
                identifier::Union{Symbol,Char}=wildcard
                ) where {N1,N2}
    Subscript{N}() where N
    Subscript(opattern::NTuple{N,Int}) where N

The subscripts of some orbital/spin degrees of freedom.
"""
struct Subscript{N1,N2,I<:Tuple,O<:Tuple,M<:Union{Function,Nothing},C<:Union{Function,Nothing},D<:Union{Symbol,Char}}
    ipattern::I
    opattern::O
    mapping::M
    constrain::C
    identifier::D
    function Subscript( ipattern::NTuple{N1,Any},
                        opattern::NTuple{N2,Any},
                        mapping::Union{Function,Nothing}=nothing,
                        constrain::Union{Function,Nothing}=nothing,
                        identifier::Union{Symbol,Char}=wildcard
                        ) where {N1,N2}
        new{N1,N2,typeof(ipattern),typeof(opattern),typeof(mapping),typeof(constrain),typeof(identifier)}(ipattern,opattern,mapping,constrain,identifier)
    end
end
Subscript{N}() where N=Subscript((wildcard,),NTuple{N,Char}(wildcard for i=1:N),nothing,nothing,wildcard)
Subscript(opattern::NTuple{N,Int}) where N=Subscript((),opattern,nothing,nothing,constant)

"""
    ==(sub1::Subscript,sub2::Subscript) -> Bool
    isequal(sub1::Subscript,sub2::Subscript) -> Bool

Judge whether two subscripts are equivalent to each other.
"""
Base.:(==)(sub1::Subscript,sub2::Subscript)=sub1.ipattern==sub2.ipattern && sub1.opattern==sub2.opattern && sub1.identifier==sub2.identifier
Base.isequal(sub1::Subscript,sub2::Subscript)=isequal(sub1.ipattern,sub2.ipattern) && isequal(sub1.opattern,sub2.opattern) && isequal(sub1.identifier,sub2.identifier)

"""
    show(io::IO,subscript::Subscript)

Show a subscript.
"""
function Base.show(io::IO,subscript::Subscript)
    if subscript.identifier==constant
        @printf io "(%s)" join(subscript.opattern,',')
    elseif subscript.identifier==wildcard
        @printf io "%s=>(%s)" join(subscript.ipattern,',') join(subscript.opattern,',')
    else
        @printf io "(%s)=>(%s) with %s" join(subscript.ipattern,',') join(subscript.opattern,',') subscript.identifier
    end
end

"""
    rank(subscript::Subscript) -> Int
    rank(::Type{<:Subscript{N}}) where N -> Int

Get the number of the independent variables that are used to describe the subscripts of some orbital/spin degrees of freedom.
"""
rank(subscript::Subscript)=subscript|>typeof|>rank
rank(::Type{<:Subscript{N}}) where N=N

"""
    dimension(subscript::Subscript) -> Int
    dimension(::Type{<:Subscript{N1,N2}}) where {N1,N2} -> Int

Get the number of the whole variables that are used to describe the subscripts of some orbital/spin degrees of freedom.
"""
dimension(subscript::Subscript)=subscript|>typeof|>dimension
dimension(::Type{<:Subscript{N1,N2}}) where {N1,N2}=N2

"""
    (subscript::Subscript{N})(::Val{'M'},values::Vararg{Int,N}) where N -> NTuple{dimension(subscript),Int}
    (subscript::Subscript{N})(::Val{'C'},values::Vararg{Int,N}) where N -> Bool

* Construct the subscripts from a set of independent variables.
* Judge whether a set of independent variables are valid to construct the subscripts.
"""
function (subscript::Subscript)() end
(subscript::Subscript{N})(::Val{'M'},values::Vararg{Int,N}) where N=subscript.mapping(values...)
(subscript::Subscript{0,N,<:Tuple,<:Tuple,Nothing,Nothing})(::Val{'M'}) where N=subscript.opattern
(subscript::Subscript{1,N,<:Tuple,<:Tuple,Nothing,Nothing})(::Val{'M'},value::Int) where N=NTuple{N,Int}(value for i=1:N)
(subscript::Subscript{N})(::Val{'C'},values::Vararg{Int,N}) where N=subscript.constrain(values...)
(subscript::Subscript{N1,N2,<:Tuple,<:Tuple,<:Union{Function,Nothing},Nothing})(::Val{'C'},values::Vararg{Int,N1}) where {N1,N2}=true

"""
    @subscript expr::Expr with constrain::Expr -> Subscript

Construct a subscript from a map and optionally with a constrain.
"""
macro subscript(expr::Expr,with::Symbol=:with,constrain::Union{Expr,Symbol}=:nothing,gensym::Expr=:(gensym=false))
    @assert expr.head==:call && expr.args[1]==:(=>)
    @assert isa(expr.args[2],Expr) && expr.args[2].head==:tuple && all(expr.args[2].args.≠Symbol(wildcard)) && all(expr.args[2].args.≠Symbol(constant))
    @assert isa(expr.args[3],Expr) && expr.args[3].head==:tuple && all(expr.args[3].args.≠Symbol(wildcard)) && all(expr.args[3].args.≠Symbol(constant))
    @assert with==:with
    @assert isa(gensym,Expr) && gensym.head==:(=) && gensym.args[1]==:gensym && isa(gensym.args[2],Bool)
    ip=gensym.args[2] ? Tuple(Base.gensym(arg) for arg in expr.args[2].args) : Tuple(expr.args[2].args)
    op=gensym.args[2] ? Tuple(ip[findfirst(isequal(arg),expr.args[2].args)] for arg in expr.args[3].args) : Tuple(expr.args[3].args)
    mapname,identifier=Base.gensym(),QuoteNode(Base.gensym())
    if constrain===:nothing
        return quote
            $mapname($(expr.args[2].args...))=$(expr.args[3])
            Subscript($ip,$op,$mapname,nothing,$identifier)
        end
    else
        constrainname=Base.gensym()
        return quote
            $mapname($(expr.args[2].args...))=$(expr.args[3])
            $constrainname($(expr.args[2].args...))=$(constrain)
            Subscript($ip,$op,$mapname,$constrainname,$identifier)
        end
    end
end

"""
    Subscripts(contents::Subscript...)

A complete set of all the independent subscripts of the orbital/spin degrees of freedom.
"""
struct Subscripts{T<:Tuple,R,D} <: CompositeTuple{T}
    contents::T
    function Subscripts(contents::T) where T<:Tuple
        R=sum(rank(fieldtype(T,i)) for i=1:fieldcount(T))
        D=sum(dimension(fieldtype(T,i)) for i=1:fieldcount(T))
        new{T,R,D}(contents)
    end
end
Subscripts(contents::Subscript...)=Subscripts(contents)

"""
    rank(subscripts::Subscripts) -> Int
    rank(::Type{S}) where S<:Subscripts -> Int

Get the total number of the independent variables of the complete subscript set.
"""
rank(subscripts::Subscripts)=subscripts|>typeof|>rank
rank(::Type{<:Subscripts{<:Tuple,R}}) where R=R

"""
    rank(subscripts::Subscripts,i::Int) -> Int
    rank(::Type{<:Subscripts{T}},i::Int) where T -> Int

Get the number of the independent variables of a component of the complete subscript set.
"""
rank(subscripts::Subscripts,i::Int)=rank(typeof(subscripts),i)
rank(::Type{<:Subscripts{T}},i::Int) where T=rank(fieldtype(T,i))

"""
    dimension(subscripts::Subscripts) -> Int
    dimension(::Type{S}) where S<:Subscripts -> Int

Get the total number of the whole variables of the complete subscript set.
"""
dimension(subscripts::Subscripts)=subscripts|>typeof|>dimension
dimension(::Type{<:Subscripts{<:Tuple,R,D}}) where {R,D}=D

"""
    dimension(subscripts::Subscripts,i::Int) -> Int
    dimension(::Type{<:Subscripts{T}},i::Int) where T -> Int

Get the total number of the whole variables of a component of the complete subscript set.
"""
dimension(subscripts::Subscripts,i::Int)=dimension(typeof(subscripts),i)
dimension(::Type{<:Subscripts{T}},i::Int) where T=dimension(fieldtype(T,i))

"""
    (subscripts::Subscripts)(::Val{'M'},values::NTuple{N,Int}) where N -> NTuple{dimension(subscripts),Int}
    (subscripts::Subscripts)(::Val{'C'},values::NTuple{N,Int}) where N -> Bool

* Construct the complete set of subscripts from a complete set of independent variables.
* Judge whether a complete set of independent variables are valid to construct the complete subscripts.
"""
function (subscripts::Subscripts)() end
@generated function (subscripts::Subscripts)(::Val{'M'},values::NTuple{N,Int}) where N
    @assert rank(subscripts)==N
    exprs,count=[],1
    for i=1:length(subscripts)
        vs=Tuple(:(values[$j]) for j=count:(count+rank(subscripts,i)-1))
        push!(exprs,:(subscripts.contents[$i](Val('M'),$(vs...))...))
        count=count+rank(subscripts,i)
    end
    return Expr(:tuple,exprs...)
end
@generated function (subscripts::Subscripts)(::Val{'C'},values::NTuple{N,Int}) where N
    @assert rank(subscripts)==N
    length(subscripts)==0 && return :(true)
    count=rank(subscripts,1)
    vs=Tuple(:(values[$j]) for j=1:count)
    expr=:(subscripts.contents[1](Val('C'),$(vs...)))
    for i=2:length(subscripts)
        vs=Tuple(:(values[$j]) for j=(count+1):(count+rank(subscripts,i)))
        expr=Expr(:(&&),expr,:(subscripts.contents[$i](Val('C'),$(vs...))))
        count=count+rank(subscripts,i)
    end
    return expr
end

"""
    *(sub1::Subscript,sub2::Subscript) -> Subscripts
    *(subs::Subscripts,sub::Subscript) -> Subscripts
    *(sub::Subscript,subs::Subscripts) -> Subscripts
    *(subs1::Subscripts,subs2::Subscripts) -> Subscripts

Get the multiplication between subscripts or complete sets of subscripts.
"""
Base.:*(sub1::Subscript,sub2::Subscript)=Subscripts(sub1,sub2)
Base.:*(subs::Subscripts,sub::Subscript)=Subscripts(subs.contents...,sub)
Base.:*(sub::Subscript,subs::Subscripts)=Subscripts(sub,subs.contents...)
Base.:*(subs1::Subscripts,subs2::Subscripts)=Subscripts(subs1.contents...,subs2.contents...)

"""
    expand(subscripts::Subscripts,dimensions::NTuple{N,Int}) where N -> SbExpand

Expand a complete set of subscripts with a given set of variable ranges.
"""
function expand(subscripts::Subscripts,dimensions::NTuple{N,Int}) where N
    @assert dimension(subscripts)==N "expand error: dismatched input dimensions $dimensions."
    dims,dcount,rcount=Vector{Int}(undef,rank(subscripts)),0,0
    for i=1:length(subscripts)
        ipattern=subscripts[i].ipattern
        opattern=subscripts[i].opattern
        for j=1:dimension(subscripts,i)
            isa(opattern[j],Int) && @assert 0<opattern[j]<=dimensions[dcount+j] "expand error: opattern($opattern) out of range."
        end
        for j=1:rank(subscripts,i)
            index,flag=1,false
            while (index=findnext(isequal(ipattern[j]),opattern,index))≠nothing
                flag || (dims[rcount+j]=dimensions[dcount+index];flag=true)
                flag && @assert dimensions[dcount+index]==dims[rcount+j] "expand error: dismatched input dimensions."
                index=index+1
            end
        end
        dcount=dcount+dimension(subscripts,i)
        rcount=rcount+rank(subscripts,i)
    end
    return SbExpand(subscripts,NTuple{rank(subscripts),Int}(dims))
end
struct SbExpand{N,S<:Subscripts,D}
    subscripts::S
    dims::NTuple{D,Int}
    function SbExpand(subscripts::Subscripts,dims::NTuple{D,Int}) where D
        @assert D==rank(subscripts) "SbExpand error: dismatched inputs."
        new{dimension(subscripts),typeof(subscripts),D}(subscripts,dims)
    end
end
Base.eltype(::Type{<:SbExpand{N}}) where N=NTuple{N,Int}
Base.IteratorSize(::Type{<:SbExpand})=Base.SizeUnknown()
function Base.iterate(sbe::SbExpand,state::Int=1)
    while state<=prod(sbe.dims)
        inds=indtosub(sbe.dims,state,corder)
        sbe.subscripts(Val('C'),inds) && return (sbe.subscripts(Val('M'),inds),state+1)
        state=state+1
    end
    return nothing
end

"""
    Coupling{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I}

The coupling intra/inter interanl degrees of freedom at different lattice points.
"""
abstract type Coupling{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I} end

defaultcenter(::Type{<:Coupling},i::Int,n::Int,::Val{1})=1
defaultcenter(::Type{<:Coupling},i::Int,n::Int,::Val{R}) where R=error("defaultcenter error: no default center for a rank-$R bond.")
@generated function propercenters(::Type{C},centers::NTuple{N,Any},::Val{R}) where {C<:Coupling,N,R}
    exprs=[:(isa(centers[$i],Int) ? (0<centers[$i]<=R ? centers[$i] : error("propercenters error: center out of range.")) : defaultcenter(C,$i,N,Val(R))) for i=1:N]
    return Expr(:tuple,exprs...)
end

"""
    Couplings{I<:ID,C<:Coupling} <: AbstractDict{I,C}

A pack of couplings intra/inter interanl degrees of freedom at different lattice points.

Alias for `Elements{I,C}`.
"""
const Couplings{I<:ID,C<:Coupling}=Elements{I,C}
Couplings(cps::Coupling...)=Elements(cps...)

end #module
