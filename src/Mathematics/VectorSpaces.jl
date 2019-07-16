module VectorSpaces

using ..Combinatorics: AbstractCombinatorics
using ...Prerequisites: rawtype
using ...Prerequisites.TypeTraits: efficientoperations,forder,corder,subtoind,indtosub

import ...Interfaces: dimension,⊕,rank,dims,inds

export VectorSpace
export HasTable,TableSorted,IsMultiIndexable,MultiIndexOrderStyle
export SimpleVectorSpace
export OrderedIndices,SimpleIndices,TabledIndices
export NamedVectorSpace

"""
    VectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type VectorSpace{B} <: AbstractVector{B} end
Base.:(==)(vs1::VectorSpace,vs2::VectorSpace) = ==(efficientoperations,vs1,vs2)
Base.isequal(vs1::VectorSpace,vs2::VectorSpace)=isequal(efficientoperations,vs1,vs2)
Base.size(vs::VectorSpace)=(vs|>dimension,)
Base.iterate(vs::VectorSpace,state::Integer=1)=state>length(vs) ? nothing : (vs[state],state+1)

"""
    HasTable(B::Bool)
    HasTable(::Type{<:VectorSpace})

Trait of whether a subtype of `VectorSpace` has the attribute `:table`.

Only two instances are allowed, the first of which is the default for a subtype:
* `HasTable(false)`: indication of not having the attribute `:table`
* `HasTable(true)`: indication of having the attribute `:table`
"""
struct HasTable{B} end
HasTable(B::Bool)=HasTable{B}()
HasTable(::Type{<:VectorSpace})=HasTable{false}()

"""
    TableSorted(B::Bool)
    TableSorted(::Type{<:VectorSpace})

Trait of whether the attribute `:table` of a subtype of `VectorSpace` is sorted.

Only two instances are allowed, the first of which is the default for a subtype:
* `TableSorted(false)`: indication of unsorted attribute `:table`
* `TableSorted(true)`: indication of sorted attribute `:table`
"""
struct TableSorted{B} end
TableSorted(B::Bool)=TableSorted{B}()
TableSorted(::Type{<:VectorSpace})=TableSorted{false}()

"""
    IsMultiIndexable(B::Bool)
    IsMultiIndexable(::Type{<:VectorSpace})

Trait of whether the bases of a subtype of `VectorSpace` can be represented by multiindices (Cartesian indices).

Only two instances are allowed, the first of which is the default for a subtype:
* `IsMultiIndexable(false)`: indication of irrepresentability by Cartesian indices
* `IsMultiIndexable(true)`: indication of representability by Cartesian indices
"""
struct IsMultiIndexable{B} end
IsMultiIndexable(B::Bool)=IsMultiIndexable{B}()
IsMultiIndexable(::Type{<:VectorSpace})=IsMultiIndexable{false}()

"""
    MultiIndexOrderStyle(M::Char)
    MultiIndexOrderStyle(::Type{<:VectorSpace})

Trait of the order style of the Cartesian-index representation of the bases of a subtype of `VectorSpace`.

Only two instances are allowed, the first of which is the default for a subtype:
* `MultiIndexOrderStyle(`C`)`: indication of C/C++ order style
* `MultiIndexOrderStyle(`F`)`: indication of Fortran order style
"""
struct MultiIndexOrderStyle{M} end
MultiIndexOrderStyle(M::Char)=M∈('C','F') ? MultiIndexOrderStyle{M}() : error("MultiIndexOrderStyle error: wrong input M($M).")
MultiIndexOrderStyle(::Type{<:VectorSpace})=MultiIndexOrderStyle{'C'}()

"""
    dimension(vs::VectorSpace) -> Int

The dimension of a vector space.
"""
dimension(vs::VectorSpace)=_dimension_(vs,vs|>typeof|>HasTable)
_dimension_(vs::VectorSpace,::HasTable{true})=length(vs.table)
_dimension_(vs::VectorSpace,::HasTable{false})=_dimension_(vs,vs|>typeof|>IsMultiIndexable)
_dimension_(vs::VectorSpace,::IsMultiIndexable{true})=prod(dims(vs))

"""
    getindex(vs::VectorSpace,i::Int)

Get the ith basis of a vector space by the `[]` operator.
"""
Base.getindex(vs::VectorSpace,i::Int)=_getindex_(vs,i,vs|>typeof|>HasTable,vs|>typeof|>IsMultiIndexable)
_getindex_(vs::VectorSpace,i::Int,::HasTable{true},::IsMultiIndexable{false})=vs.table[i]
_getindex_(vs::VectorSpace,i::Int,::HasTable{true},::IsMultiIndexable{true})=rawtype(eltype(vs))(vs.table[i],vs)
_getindex_(vs::VectorSpace,i::Int,::HasTable{false},::IsMultiIndexable{true})=_getindex_(vs,i,vs|>typeof|>MultiIndexOrderStyle)
_getindex_(vs::VectorSpace,i::Int,::MultiIndexOrderStyle{'F'})=rawtype(eltype(vs))(indtosub(dims(vs),i,forder),vs)
_getindex_(vs::VectorSpace,i::Int,::MultiIndexOrderStyle{'C'})=rawtype(eltype(vs))(indtosub(dims(vs),i,corder),vs)

"""
    searchsortedfirst(vs::VectorSpace{B},basis::B) where B -> Int
    searchsortedfirst(vs,basis) -> Int

Search the index of a basis in a vector space.
"""
Base.searchsortedfirst(vs::VectorSpace{B},basis::B) where B=_search_(vs,basis,vs|>typeof|>HasTable)
_search_(vs::VectorSpace{B},basis::B,::HasTable{true}) where B=_search_(vs,basis,vs|>typeof|>TableSorted,vs|>typeof|>IsMultiIndexable)
_search_(vs::VectorSpace{B},basis::B,::TableSorted{true},::IsMultiIndexable{false}) where B=searchsortedfirst(vs.table,basis)
_search_(vs::VectorSpace{B},basis::B,::TableSorted{false},::IsMultiIndexable{false}) where B=findfirst(isequal(basis),vs.table)
_search_(vs::VectorSpace{B},basis::B,::TableSorted{true},::IsMultiIndexable{true}) where B=searchsortedfirst(vs.table,inds(basis,vs))
_search_(vs::VectorSpace{B},basis::B,::TableSorted{false},::IsMultiIndexable{true}) where B=findfirst(isequal(inds(basis,vs)),vs.table)
_search_(vs::VectorSpace{B},basis::B,::HasTable{false}) where B=_search_(vs,basis,vs|>typeof|>IsMultiIndexable)
_search_(vs::VectorSpace{B},basis::B,::IsMultiIndexable{true}) where B=_search_(vs,basis,vs|>typeof|>MultiIndexOrderStyle)
function _search_(vs::VectorSpace{B},basis::B,::MultiIndexOrderStyle{'F'}) where B
    basis ∈ vs ? subtoind(dims(vs),inds(basis,vs),forder) : error("searchsortedfirst error: input basis($basis) not found.")
end
function _search_(vs::VectorSpace{B},basis::B,::MultiIndexOrderStyle{'C'}) where B
    basis ∈ vs ? subtoind(dims(vs),inds(basis,vs),corder) : error("searchsortedfirst error: input basis($basis) not found.")
end
function _search_(vs::VectorSpace{B},basis::B,::IsMultiIndexable{false}) where B
    for i=1:dimension(vs)
        vs[i]==basis && return i
    end
end
function Base.searchsortedfirst(vs,basis)
    lo,hi=0,length(vs)+1
    @inbounds while lo<hi-1
        m=(lo+hi)>>>1
        if vs[m]<basis
            lo=m
        else
            hi=m
        end
    end
    return hi
end

"""
    findfirst(basis::B,vs::VectorSpace{B}) where B -> Int
    findfirst(bases,vs::VectorSpace) -> NTuple{length(bases),Int}

Get the index of a basis or the indices of a couple of bases in a vector space.
"""
Base.findfirst(basis::B,vs::VectorSpace{B}) where B=_find_(basis,vs,vs|>typeof|>HasTable)
Base.findfirst(bases,vs::VectorSpace)=NTuple{length(bases),Int}(findfirst(basis,vs)::Int for basis in bases)
_find_(basis::B,vs::VectorSpace{B},::HasTable{true}) where B=_find_(basis,vs,vs|>typeof|>TableSorted)
_find_(basis::B,vs::VectorSpace{B},::TableSorted{true}) where B=_findwithmissing_(basis,vs)
_find_(basis::B,vs::VectorSpace{B},::TableSorted{false}) where B=_findwithoutmissing_(basis,vs)
_find_(basis::B,vs::VectorSpace{B},::HasTable{false}) where B=_find_(basis,vs,vs|>typeof|>IsMultiIndexable)
_find_(basis::B,vs::VectorSpace{B},::IsMultiIndexable{true}) where B=_findwithoutmissing_(basis,vs)
_find_(basis::B,vs::VectorSpace{B},::IsMultiIndexable{false}) where B=_findwithmissing_(basis,vs)
_findwithmissing_(basis::B,vs::VectorSpace{B}) where B=(index=searchsortedfirst(vs,basis);0<index<=dimension(vs) && vs[index]==basis ? index : nothing)
_findwithoutmissing_(basis::B,vs::VectorSpace{B}) where B=searchsortedfirst(vs,basis)

"""
    in(basis::B,vs::VectorSpace{B}) where B -> Bool

Judge whether a basis belongs to a vector space.
"""
Base.in(basis::B,vs::VectorSpace{B}) where B=_in_(basis,vs,vs|>typeof|>HasTable,vs|>typeof|>IsMultiIndexable)
_in_(basis::B,vs::VectorSpace{B},::HasTable,::IsMultiIndexable) where B=isa(findfirst(basis,vs),Integer)
_in_(basis::B,vs::VectorSpace{B},::HasTable{false},::IsMultiIndexable{true}) where B=_in_(basis,vs,Val(vs|>typeof|>rank))
@generated function _in_(basis::B,vs::VectorSpace{B},::Val{N}) where {B,N}
    N==0 && return :(true)
    prepare=:(INDS=inds(basis,vs);DIMS=dims(vs))
    expr=:(0<INDS[1]<=DIMS[1])
    for i=2:N
        expr=Expr(:&&,expr,:(0<INDS[$i]<=DIMS[$i]))
    end
    return Expr(:block,prepare,expr)
end

"""
    SimpleVectorSpace{S}(table::NTuple{N,B}) where {S,B,N}
    SimpleVectorSpace{S}(table...) where S

Simplest vector space, whose bases are stored in the attribute `:table` as an ntuple.

The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted.
"""
struct SimpleVectorSpace{S,B,N} <: VectorSpace{B}
    table::NTuple{N,B}
    SimpleVectorSpace{S}(table::NTuple{N,B}) where {S,B,N}=(@assert(S∈('T','F'),"SimpleVectorSpace error: wrong type parameter S($S).");new{S,B,N}(table))
end
SimpleVectorSpace{S}(table...) where S=SimpleVectorSpace{S}(table)
HasTable(::Type{<:SimpleVectorSpace})=HasTable(true)
TableSorted(::Type{<:SimpleVectorSpace{'T'}})=TableSorted(true)
TableSorted(::Type{<:SimpleVectorSpace{'F'}})=TableSorted(false)

"""
    ⊕(basis1::B,basis2::B) where B -> SimpleVectorSpace{'F',B,2}
    ⊕(basis::B,vs::SimpleVectorSpace{<:Any,B}) where B -> SimpleVectorSpace{'F',B}
    ⊕(vs::SimpleVectorSpace{<:Any,B},basis::B) where B -> SimpleVectorSpace{'F',B}
    ⊕(vs1::SimpleVectorSpace{<:Any,B},vs2::SimpleVectorSpace{<:Any,B}) where B -> SimpleVectorSpace{'F',B}

Get the direct sum between bases or simple vector spaces.
"""
⊕(basis1::B,basis2::B) where B=SimpleVectorSpace{'F'}(basis1,basis2)
⊕(basis::B,vs::SimpleVectorSpace{<:Any,B}) where B=SimpleVectorSpace{'F'}(basis,vs.table...)
⊕(vs::SimpleVectorSpace{<:Any,B},basis::B) where B=SimpleVectorSpace{'F'}(vs.table...,basis)
⊕(vs1::SimpleVectorSpace{<:Any,B},vs2::SimpleVectorSpace{<:Any,B}) where B=SimpleVectorSpace{'F'}(vs1.table...,vs2.table...)

"""
    OrderedIndices{N} <: VectorSpace{NTuple{N,Int}}

The simplest abstract class of multiindexable vector spaces, whose bases are just tuples of integers.

This class of vector spaces must have the following attribute:
`dims::NTuple{N,Int}`: the dimesnions of the Cartesian indices along all axes
"""
abstract type OrderedIndices{N} <: VectorSpace{NTuple{N,Int}} end
IsMultiIndexable(::Type{<:OrderedIndices})=IsMultiIndexable(true)
dims(vs::OrderedIndices)=vs.dims
inds(basis::NTuple{N,Int},::OrderedIndices{N}) where N=basis
rank(::Type{<:OrderedIndices{N}}) where N=N
Tuple(index::NTuple{N,Int},::OrderedIndices{N}) where N=index

"""
    SimpleIndices{M}(dims::NTuple{N,Int}) where {M,N}

Simple ordered Cartesian indices.

!!! note
    1) It can be C/C++ ordered or Fortran ordered depending on the first type parameter `M`, with `'C'` for the former and `'F'` the latter.
    2) For its bases (Cartesian indices), there is no restriction except that they should be in the proper range defined by its `dims`.
"""
struct SimpleIndices{M,N} <: OrderedIndices{N}
    dims::NTuple{N,Int}
    SimpleIndices{M}(dims::NTuple{N,Int}) where {M,N}=(@assert(M∈('C','F'),"SimpleIndices error: wrong type parameter M($M).");new{M,N}(dims))
end
SimpleIndices{M}(dims::Int...) where M=SimpleIndices{M}(dims)
SimpleIndices{M,N}(dims::NTuple{N,Int}) where {M,N}=SimpleIndices{M}(dims)
MultiIndexOrderStyle(::Type{<:SimpleIndices{M}}) where M=MultiIndexOrderStyle(M)

"""
    TabledIndices{S}(dims::NTuple{N,Int},table::Vector{NTuple{N,Int}}) where {S,N}
    TabledIndices{N}(::Type{M},len::Int) where {N,M<:AbstractCombinatorics}

Tabled ordered Cartesian indices.

Compared to [`SimpleIndices`](@ref), the bases of this kind of vector spaces are stored in the attribute `:table`, which must be a vector of tuple of integers. The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted. This type suits the situations when the Cartesian indices are restricted by extra conditions except that propoesed by the attribute `:dims`.
"""
struct TabledIndices{S,N} <: OrderedIndices{N}
    dims::NTuple{N,Int}
    table::Vector{NTuple{N,Int}}
    TabledIndices{S}(dims::NTuple{N,Int},table::Vector{NTuple{N,Int}}) where {S,N}=new{S,N}(dims,table)
end
TabledIndices{N}(::Type{M},len::Int) where {N,M<:AbstractCombinatorics}=TabledIndices{'T'}(ntuple(i->len,N),M{N}(1:len)|>collect)
HasTable(::Type{<:TabledIndices})=HasTable(true)
TableSorted(::Type{<:TabledIndices{'T'}})=TableSorted(true)
TableSorted(::Type{<:TabledIndices{'F'}})=TableSorted(false)

"""
    NamedVectorSpace{M,NS,BS<:Tuple,VS<:Tuple{Vararg{AbstractVector}}} <: VectorSpace{NamedTuple{NS,BS}}

Abstract named vector space.

This is a wrapper of multiindexable vector spaces, each of whose indexable dimensions is associated with a name.

It has four type parameters:
* `M`: mode of the named vector space. It specifies how the indexable dimensions are combined to form the bases of the named vector space, and must take one of the following values:
  - `:zip`: elements from each indexable dimensions are zipped together to form the bases,
  - `:product`: elements from each indexable dimensions are direct producted together to form the bases.
For the `:zip` mode, all the indexable dimensions should have the same number of elements, and the number of formed bases is equal to this number; for the `:product` mode, there are no restriction on the numbers of the indexable dimensions, and the number of the final bases is equal to their product.
* `NS::Tuple{Vararg{Symbol}}`: the names of the indexable dimensions
* `BS<:Tuple`: the eltypes of the indexable dimensions
* `VS<:Tuple{Vararg{AbstractVector}}`: the contents of the indexable dimensions

The concrete types must have the following attribute:
* `:contents::VS`: storage of the contents of the indexable dimensions

By default, a named vector space uses C order for the indexable dimensions when the mode is `:product`. You can change it to F order for your own subtypes by defining the [`MultiIndexOrderStyle`](@ref) trait.
"""
abstract type NamedVectorSpace{M,NS,BS<:Tuple,VS<:Tuple{Vararg{AbstractVector}}} <: VectorSpace{NamedTuple{NS,BS}} end
IsMultiIndexable(::Type{<:NamedVectorSpace})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{<:NamedVectorSpace})=MultiIndexOrderStyle('C')
rank(::Type{<:NamedVectorSpace{:zip}})=1
rank(::Type{<:NamedVectorSpace{:product,NS}}) where NS=length(NS)
dims(nvs::NamedVectorSpace{:zip})=(length(getfield(nvs,:contents)[1]),)
@generated dims(nvs::NamedVectorSpace{:product})=Expr(:tuple,[:(length(getfield(nvs,:contents)[$i])) for i=1:rank(nvs)]...)
function inds(basis::NamedTuple{NS},nvs::NamedVectorSpace{:zip,NS}) where NS
    index=findfirst(isequal(basis[1]),getfield(nvs,:contents)[1])
    @assert isa(index,Integer) "inds error: input basis out of range."
    for i=2:length(NS)
        @assert getfield(nvs,:contents)[i][index]==basis[i] "inds error: input basis out of range."
    end
    return (index,)
end
@generated function inds(basis::NamedTuple{NS,BS},nvs::NamedVectorSpace{:product,NS,BS}) where {NS,BS<:Tuple}
    return Expr(:tuple,[:(findfirst(isequal(basis[$i]),getfield(nvs,:contents)[$i])) for i=1:length(NS)]...)
end
@generated function NamedTuple(index::Tuple{Int},nvs::NamedVectorSpace{:zip,NS}) where NS
    return Expr(:tuple,[:($name=getfield(nvs,:contents)[$i][index[1]]) for (i,name) in enumerate(NS)]...)
end
@generated function NamedTuple(index::NTuple{N,Int},nvs::NamedVectorSpace{:product,NS}) where {N,NS}
    @assert N==length(NS) "NamedTuple error: dismatched input index and rank of named vector space."
    return Expr(:tuple,[:($name=getfield(nvs,:contents)[$i][index[$i]]) for (i,name) in enumerate(NS)]...)
end

"""
    keys(nvs::NamedVectorSpace) -> Tuple{Vararg{Symbol}}
    keys(::Type{<:NamedVectorSpace{M,NS} where M}) where NS -> Tuple{Vararg{Symbol}}

Get the names of a named vector space.
"""
Base.keys(nvs::NamedVectorSpace)=nvs|>typeof|>keys
Base.keys(::Type{<:NamedVectorSpace{M,NS} where M}) where NS=NS

"""
    values(nvs::NamedVectorSpace) -> Tuple{Vararg{AbstractVector}}

Get the contents of a named vector space.
"""
Base.values(nvs::NamedVectorSpace)=getfield(nvs,:contents)

"""
    pairs(nvs::NamedVectorSpace) -> Base.Iterators.Pairs

Get the name-value pairs of a named vector space.
"""
Base.pairs(nvs::NamedVectorSpace)=Base.Generator(=>,nvs|>keys,nvs|>values)

"""
    eltype(nvs::NamedVectorSpace,i::Int)
    eltype(::Type{<:NamedVectorSpace{M,NS,BS} where {M,NS}},i::Int) where BS

Get the eltype of the ith indexable dimensions of a named vector space.
"""
Base.eltype(nvs::NamedVectorSpace,i::Int)=eltype(nvs|>typeof,i)
Base.eltype(::Type{<:NamedVectorSpace{M,NS,BS} where {M,NS}},i::Int) where BS=fieldtype(BS,i)

end # module
