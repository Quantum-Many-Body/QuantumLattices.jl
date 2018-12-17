module DegreeOfFreedom

using Printf: @printf
using ...Utilities: Float
using ...Utilities.NamedVector: AbstractNamedVector
using ...Utilities.CompositeStructure: CompositeDict
using ...Utilities.AlgebraOverField: SimpleID, Element, Elements
using ..Spatial: PID

export IID,Index,pidtype,pid,iidtype,iid
export IndexToTuple,DirectIndexToTuple,directindextotuple,FilteredAttributes
export Internal,IDFConfig,Table,Coupling,Couplings

"""
    IID

The id of an internal degree of freedom.
"""
abstract type IID <: AbstractNamedVector end

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
    return :(pidtype(index)($(exprs...)))
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
    return :(iidtype(index)($(exprs...)))
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
(indextotuple::DirectIndexToTuple)(index::Index)=invoke(convert,Tuple{Type{Tuple},AbstractNamedVector},Tuple,index)

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
    Internal

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} end

"""
    eltype(internal::Internal)
    eltype(::Type{<:Internal{I}}) where I

Get the type of the `IID`s that an `Internal` contains.
"""
Base.eltype(internal::Internal)=internal|>typeof|>eltype
Base.eltype(::Type{<:Internal{I}}) where I=I

"""
    ==(i1::I,i2::I) where I<:Internal
    isequal(i1::I,i2::I) where I<:Internal

Compare two internals and judge whether they are equal to each other.
"""
Base.:(==)(i1::I,i2::I) where I<:Internal=all(getfield(i1,i)==getfield(i2,i) for i=1:fieldcount(I))
Base.isequal(i1::I,i2::I) where I<:Internal=all(isequal(getfield(i1,i),getfield(i2,i)) for i=1:fieldcount(I))

"""
    show(io::IO,i::Internal)

Show an internal.
"""
Base.show(io::IO,i::Internal)=@printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i,name))" for name in i|>typeof|>fieldnames),",")

"""
    IDFConfig(map::Function,::Type{I},pids::AbstractVector{<:PID}=[]) where I<:Internal

Configuration of the internal degrees of freedom at a lattice.

`map` maps a `PID` to an `Internal`.
"""
struct IDFConfig{P<:PID,I<:Internal} <: CompositeDict{P,I}
    map::Function
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

"""
    Coupling{V,I}

The coupling intra/inter interanl degrees of freedom at different lattice points.
"""
abstract type Coupling{V,I} <: Element{V,I} end

"""
    Couplings{I,C<:Coupling} <: AbstractDict{I,C}

A pack of couplings intra/inter interanl degrees of freedom at different lattice points.

Alias for `Elements{I,C}`.
"""
const Couplings{I,C<:Coupling}=Elements{I,C}

end #module
