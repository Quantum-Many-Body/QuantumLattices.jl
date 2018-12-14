module DegreeOfFreedom

using Printf: @printf
using ...Utilities: Float
using ...Utilities.NamedVector: AbstractNamedVector
using ...Utilities.CompositeStructure: CompositeDict
using ...Utilities.AlgebraOverField: SimpleID, Element, Elements
using ..Spatial: PID

export IID,Index,pidtype,pid,iidtype,iid
export IndexToTuple,DirectIndexToTuple,directindextotuple
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

The type of the spatial part of an index.
"""
pidtype(index::Index)=index|>typeof|>pidtype
pidtype(::Type{<:Index{P,I}}) where {P,I}=P

"""
    pid(index::Index) -> PID

The spatial part of an index.
"""
pid(index::Index)=pidtype(index)((getfield(index,name) for name in index|>pidtype|>fieldnames)...)

"""
    iidtype(index::Index)
    iidtype(::Type{<:Index{P,I}}) where {P,I}

The type of the internal part of an index.
"""
iidtype(index::Index)=index|>typeof|>iidtype
iidtype(::Type{<:Index{P,I}}) where {P,I}=I

"""
    iid(index::Index) -> IID

The internal part of an index.
"""
iid(index::Index)=iidtype(index)((getfield(index,name) for name in index|>iidtype|>fieldnames)...)

"""
    union(::Type{P},::Type{I}) where {P<:PID,I<:IID}

Combine a concrete `PID` type and a concrete `IID` type to a concrete `Index` type.
"""
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:IID}=Index{P,I}

"""
    IndexToTuple

The rules for converting an index to a tuple.

As a function, every instance should accept two positional arguments
* `index::Index`: the index to be converted to a tuple
* `mask::NTuple{N,Symbol}`: the names of the attributes of an index to be omitted during the conversion
"""
abstract type IndexToTuple <:Function end

"""
    DirectIndexToTuple

Direct index to tuple.
"""
struct DirectIndexToTuple <: IndexToTuple end

"""
    (indextotuple::DirectIndexToTuple)(index::Index) -> Tuple
    (indextotuple::DirectIndexToTuple)(index::Index,mask::Symbol...) -> Tuple

Convert an index to tuple directly.

When `mask` is nonempty, those attributes of `Index` whose names are in `mask` will be omitted during the conversion.
"""
(indextotuple::DirectIndexToTuple)(index::Index)=invoke(convert,Tuple{Type{Tuple},AbstractNamedVector},Tuple,index)
(indextotuple::DirectIndexToTuple)(index::Index,mask::Symbol...)=Tuple(getfield(index,name) for name in index|>typeof|>fieldnames if name ∉ mask)

"""
    directindextotuple

Indicate that the conversion from an index to a tuple is direct.
"""
const directindextotuple=DirectIndexToTuple()

"""
    convert(::Type{Tuple},index::Index;by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N -> Tuple

Convert an index to tuple.

`by` specifies the algorithm to convert an index to tuple, while `mask` contains the names of the omitted attributes of an index during the conversion.
"""
Base.convert(::Type{Tuple},index::Index;by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N=by(index,mask...)

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
    IDFConfig{I}(indextotuple::IndexToTuple,map::Function,pids::AbstractVector{<:PID}=[]) where I<:Internal

Configuration of the internal degrees of freedom at a lattice.

`map` maps a `PID` to an `Internal`, while `indextotuple` maps an `Index` to a tuple.
"""
struct IDFConfig{I<:Internal,P<:PID,T<:IndexToTuple} <: CompositeDict{P,I}
    map::Function
    indextotuple::T
    contents::Dict{P,I}
end
function IDFConfig{I}(map::Function,indextotuple::IndexToTuple,pids::AbstractVector{<:PID}=[]) where I<:Internal
    contents=Dict{pids|>eltype,I}()
    for pid in pids
        contents[pid]=map(pid)
    end
    IDFConfig(map,indextotuple,contents)
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
    Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N -> Table

Convert an sequence of indices to the corresponding index-sequence table.

The input indices will be converted to tuples by the `by` function along with the `mask` parameter, which contains the names of omitted attributes of the indices during this conversion. Then duplicates are removed and the resulting unique tuples are sorted, which determines the sequence of the input indices. Note that two indices have the same sequence if their converted tupels are equal to each other.
"""
function Table(indices::AbstractVector{<:Index};by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N
    tuples=[by(index,mask...) for index in indices]
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
    Table(config::IDFConfig;mask::NTuple{N,Symbol}=()) where N -> Table

Get the index-sequence table of the whole internal Hilbert spaces at a lattice.

When `mask` is nonempty, it contains the names of the attributes that will be omitted during the conversion from an index to tuple.
"""
function Table(config::IDFConfig;mask::NTuple{N,Symbol}=()) where N
    result=union(config|>keytype,config|>valtype|>eltype)[]
    for (pid,internal) in config
        for iid in internal
            push!(result,(result|>eltype)(pid,iid))
        end
    end
    Table(result,by=config.indextotuple,mask=mask)
end

"""
    union(tables::Table...;by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N -> Table

Unite several index-sequence tables.

See [`Table`](@ref) for more details.
"""
function Base.union(tables::Table...;by::IndexToTuple=directindextotuple,mask::NTuple{N,Symbol}=()) where N
    indices=(tables|>eltype|>keytype)[]
    for table in tables
        for index in keys(table)
            push!(indices,index)
        end
    end
    Table(indices,by=by,mask=mask)
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
    Coupling

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
