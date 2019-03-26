module DegreesOfFreedom

using Printf: @printf
using StaticArrays: SVector
using ..Spatials: PID,AbstractBond
using ...Interfaces: rank,dimension
using ...Prerequisites: Float,decimaltostr
using ...Prerequisites.CompositeStructures: CompositeDict
using ...Mathematics.VectorSpaces: VectorSpace
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements

import ..Spatials: pidtype,rcoord,icoord

export IID,Index,pid,iidtype,iid
export IndexToTuple,DirectIndexToTuple,directindextotuple,FilteredAttributes
export Internal,IDFConfig,Table
export OID,Operator,Operators,isHermitian,twist
export oidtype,otype

"""
    IID

The id of an internal degree of freedom.
"""
abstract type IID <: SimpleID end

"""
    Internal

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: VectorSpace{I} end

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
@generated function (INDEX::Type{<:Index})(pid::PID,iid::IID)
    INDEX<:UnionAll ? :(INDEX(convert(Tuple,pid)...,convert(Tuple,iid)...)) : :((INDEX.name.wrapper)(convert(Tuple,pid)...,convert(Tuple,iid)...))
end

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
    IDFConfig{I}(map::Function,pids::Union{AbstractVector{<:PID},Tuple{}}=()) where I<:Internal

Configuration of the internal degrees of freedom at a lattice.

`map` maps a `PID` to an `Internal`.
"""
struct IDFConfig{I<:Internal,M<:Function,P<:PID} <: CompositeDict{P,I}
    map::M
    contents::Dict{P,I}
end
function IDFConfig{I}(map::Function,pids::Union{AbstractVector{<:PID},Tuple{}}=()) where I<:Internal
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

@generated function oidcoord(vector::SVector{N,Float}) where N
    exprs=[:(vector[$i]===-0.0 ? 0.0 : vector[$i]) for i=1:N]
    return :(SVector($(exprs...)))
end

"""
    OID(index::Index,::Nothing,::Nothing,seq::Union{Nothing,Int})
    OID(index::Index,rcoord::SVector{N,Float},icoord::SVector{N,Float},seq::Union{Nothing,Int}) where N
    OID(index::Index,rcoord::Vector{Float},icoord::Vector{Float},seq::Union{Nothing,Int})
    OID(index::Index;rcoord::Union{Nothing,SVector,Vector{Float}}=nothing,icoord::Union{Nothing,SVector,Vector{Float}}=nothing,seq::Union{Nothing,Int}=nothing)

Operator id.
"""
struct OID{I<:Index,RC<:Union{Nothing,SVector},IC<:Union{Nothing,SVector},S<:Union{Nothing,Int}} <: SimpleID
    index::I
    rcoord::RC
    icoord::IC
    seq::S
    OID(index::Index,::Nothing,::Nothing,seq::Union{Nothing,Int})=new{typeof(index),Nothing,Nothing,typeof(seq)}(index,nothing,nothing,seq)
    function OID(index::Index,rcoord::SVector{N,Float},icoord::SVector{N,Float},seq::Union{Nothing,Int}) where N
        new{typeof(index),SVector{N,Float},SVector{N,Float},typeof(seq)}(index,oidcoord(rcoord),oidcoord(icoord),seq)
    end
end
OID(index::Index,rcoord::Vector{Float},icoord::Vector{Float},seq::Union{Nothing,Int})=OID(index,SVector{length(rcoord)}(rcoord),SVector{length(icoord)}(icoord),seq)
function OID(index::Index;rcoord::Union{Nothing,SVector,Vector{Float}}=nothing,icoord::Union{Nothing,SVector,Vector{Float}}=nothing,seq::Union{Nothing,Int}=nothing)
    OID(index,rcoord,icoord,seq)
end
totuple(::Nothing)=nothing
@generated totuple(v::SVector{N}) where N=Expr(:tuple,[:(v[$i]) for i=1:N]...)
Base.hash(oid::OID{<:Index},h::UInt)=hash((oid.index,totuple(oid.rcoord)),h)
Base.fieldnames(::Type{<:OID})=(:index,:rcoord,:icoord,:seq)
Base.propertynames(::Type{<:ID{<:NTuple{N,OID}}},private::Bool=false) where N=private ? (:contents,:indexes,:rcoords,:icoords) : (:indexes,:rcoords,:icoords)

"""
    show(io::IO,oid::OID)

Show an operator id.
"""
function Base.show(io::IO,oid::OID)
    @printf io "OID(%s" oid.index
    oid.rcoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.rcoord,","))
    oid.icoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.icoord,","))
    @printf io ",%s)" oid.seq===nothing ? ":" : oid.seq
end

"""
    adjoint(oid::OID) -> typeof(oid)
    adjoint(oid::ID{<:NTuple{N,OID}}) where N -> typeof(oid)

Get the adjoint of an operator id.
"""
Base.adjoint(oid::OID)=OID(oid.index',oid.rcoord,oid.icoord,oid.seq)
@generated Base.adjoint(oid::ID{<:NTuple{N,OID}}) where N=Expr(:call,:ID,[:(oid[$i]') for i=N:-1:1]...)

"""
    isHermitian(oid::ID{<:NTuple{N,OID}}) where N -> Bool

Judge whether an operator id is Hermitian.
"""
function isHermitian(oid::ID{<:NTuple{N,OID}}) where N
    for i=1:((N+1)÷2)
        oid[i]'==oid[N+1-i] || return false
    end
    return true
end

"""
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing})
    oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table})

Get the compatible oid type.
"""
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{Nothing})=OID{union(B|>pidtype,I),SVector{B|>dimension,Float},SVector{B|>dimension,Float},Nothing}
oidtype(I::Type{<:IID},B::Type{<:AbstractBond},::Type{<:Table})=OID{union(B|>pidtype,I),SVector{B|>dimension,Float},SVector{B|>dimension,Float},Int}

"""
    Operator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Element{N,V,I}

Abstract type for an operator.
"""
abstract type Operator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Element{N,V,I} end
function (O::Type{<:Operator})( value::Number,
                                indexes::NTuple{N,Index};
                                rcoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                icoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                seqs::Union{Nothing,NTuple{N,Int}}=nothing
                                ) where {N,M}
    rcoords===nothing && (rcoords=ntuple(i->nothing,N))
    icoords===nothing && (icoords=ntuple(i->nothing,N))
    seqs===nothing && (seqs=ntuple(i->nothing,N))
    return O(value,ID(OID,indexes,rcoords,icoords,seqs))
end

"""
    show(io::IO,opt::Operator)

Show an operator.
"""
Base.show(io::IO,opt::Operator)=@printf io "%s(value=%s,id=%s)" nameof(typeof(opt)) decimaltostr(opt.value) opt.id

"""
    adjoint(opt::Operator{N}) where N -> Operator

Get the adjoint of an operator.
"""
Base.adjoint(opt::Operator)=typeof(opt).name.wrapper(opt.value',opt.id')

"""
    isHermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
isHermitian(opt::Operator)=isa(opt.value,Real) && isHermitian(opt.id)

"""
    rcoord(opt::Operator{1}) -> SVector
    rcoord(opt::Operator{2}) -> SVector

Get the whole rcoord of an operator.
"""
rcoord(opt::Operator{1})=opt.id[1].rcoord
rcoord(opt::Operator{2})=opt.id[1].rcoord-opt.id[2].rcoord

"""
    icoord(opt::Operator{1}) -> SVector
    icoord(opt::Operator{2}) -> SVector

Get the whole icoord of an operator.
"""
icoord(opt::Operator{1})=opt.id[1].icoord
icoord(opt::Operator{2})=opt.id[1].icoord-opt.id[2].icoord

"""
    otype

Get the compatible operator type from a term type, a bond type and a table type.
"""
function otype end

"""
    Operators(opts::Operator...)

A set of operators.

Type alias of `Operators{I<:ID,O<:Operator}=Elements{I,O}`.
"""
const Operators{I<:ID,O<:Operator}=Elements{I,O}
Operators(opts::Operator...)=Elements(opts...)

"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
function Base.adjoint(opts::Operators)
    result=Operators{opts|>keytype,opts|>valtype}()
    for opt in values(opts)
        nopt=opt|>adjoint
        result[nopt.id]=nopt
    end
    return result
end

"""
    isHermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
isHermitian(opts::Operators)=opts==opts'

"""
    twist(operator::Operator,vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float}) -> Operator

Twist an operator.
"""
function twist(operator::Operator,vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    phase=one(operator.value)
    for i=1:rank(operator)
        phase=phase*twist(operator.id[i],vectors,values)
    end
    return replace(operator,value=operator.value*phase)
end

end #module
