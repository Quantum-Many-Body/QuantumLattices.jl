module VectorSpaces

using ..Combinatorics: AbstractCombinatorics
using ...Prerequisites.TypeTraits: efficientoperations,forder,corder,subtoind,indtosub

import ...Prerequisites.Interfaces: dimension,⊕,rank,dims,inds

export dimension,⊕,rank,dims,inds
export AbstractVectorSpace
export HasTable,TableSorted,IsMultiIndexable,MultiIndexOrderStyle
export DirectVectorSpace
export AbstractOrderedIndices,DirectOrderedIndices,TabledOrderedIndices

"""
    AbstractVectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type AbstractVectorSpace{B} <: AbstractVector{B} end
Base.:(==)(vs1::AbstractVectorSpace,vs2::AbstractVectorSpace) = ==(efficientoperations,vs1,vs2)
Base.isequal(vs1::AbstractVectorSpace,vs2::AbstractVectorSpace)=isequal(efficientoperations,vs1,vs2)
Base.size(vs::AbstractVectorSpace)=(vs|>dimension,)
Base.iterate(vs::AbstractVectorSpace,state::Integer=1)=state>length(vs) ? nothing : (vs[state],state+1)

"""
    HasTable(B::Bool)
    HasTable(::Type{<:AbstractVectorSpace})

Trait of whether a subtype of `AbstractVectorSpace` has the attribute `:table`.

Only two instances are allowed, the first of which is the default for a subtype:
* `HasTable(false)`: indication of not having the attribute `:table`
* `HasTable(true)`: indication of having the attribute `:table`
"""
struct HasTable{B} end
HasTable(B::Bool)=HasTable{B}()
HasTable(::Type{<:AbstractVectorSpace})=HasTable{false}()

"""
    TableSorted(B::Bool)
    TableSorted(::Type{<:AbstractVectorSpace})

Trait of whether the attribute `:table` of a subtype of `AbstractVectorSpace` is sorted.

Only two instances are allowed, the first of which is the default for a subtype:
* `TableSorted(false)`: indication of unsorted attribute `:table`
* `TableSorted(true)`: indication of sorted attribute `:table`
"""
struct TableSorted{B} end
TableSorted(B::Bool)=TableSorted{B}()
TableSorted(::Type{<:AbstractVectorSpace})=TableSorted{false}()

"""
    IsMultiIndexable(B::Bool)
    IsMultiIndexable(::Type{<:AbstractVectorSpace})

Trait of whether the bases of a subtype of `AbstractVectorSpace` can be represented by multiindices (Cartesian indices).

Only two instances are allowed, the first of which is the default for a subtype:
* `IsMultiIndexable(false)`: indication of irrepresentability by Cartesian indices
* `IsMultiIndexable(true)`: indication of representability by Cartesian indices
"""
struct IsMultiIndexable{B} end
IsMultiIndexable(B::Bool)=IsMultiIndexable{B}()
IsMultiIndexable(::Type{<:AbstractVectorSpace})=IsMultiIndexable{false}()

"""
    MultiIndexOrderStyle(M::Char)
    MultiIndexOrderStyle(::Type{<:AbstractVectorSpace})

Trait of the order style of the Cartesian-index representation of the bases of a subtype of `AbstractVectorSpace`.

Only two instances are allowed, the first of which is the default for a subtype:
* `MultiIndexOrderStyle(`C`)`: indication of C/C++ order style
* `MultiIndexOrderStyle(`F`)`: indication of Fortran order style
"""
struct MultiIndexOrderStyle{M} end
MultiIndexOrderStyle(M::Char)=M∈('C','F') ? MultiIndexOrderStyle{M}() : error("MultiIndexOrderStyle error: wrong input M($M).")
MultiIndexOrderStyle(::Type{<:AbstractVectorSpace})=MultiIndexOrderStyle{'C'}()

"""
    dimension(vs::AbstractVectorSpace)

The dimension of a vector space.
"""
dimension(vs::AbstractVectorSpace)=_dimension_(vs,vs|>typeof|>HasTable)
_dimension_(vs::AbstractVectorSpace,::HasTable{true})=length(vs.table)
_dimension_(vs::AbstractVectorSpace,::HasTable{false})=_dimension_(vs,vs|>typeof|>IsMultiIndexable)
_dimension_(vs::AbstractVectorSpace,::IsMultiIndexable{true})=prod(dims(vs))

"""
    getindex(vs::AbstractVectorSpace,i::Int)

Get the ith basis of a vector space by the `[]` operator.
"""
Base.getindex(vs::AbstractVectorSpace,i::Int)=_index_(vs,i,vs|>typeof|>HasTable,vs|>typeof|>IsMultiIndexable)
_index_(vs::AbstractVectorSpace,i::Int,::HasTable{true},::IsMultiIndexable{false})=vs.table[i]
_index_(vs::AbstractVectorSpace,i::Int,::HasTable{true},::IsMultiIndexable{true})=eltype(vs).name.wrapper(vs.table[i],vs)
_index_(vs::AbstractVectorSpace,i::Int,::HasTable{false},::IsMultiIndexable{true})=_index_(vs,i,vs|>typeof|>MultiIndexOrderStyle)
_index_(vs::AbstractVectorSpace,i::Int,::MultiIndexOrderStyle{'F'})=eltype(vs).name.wrapper(indtosub(dims(vs),i,forder),vs)
_index_(vs::AbstractVectorSpace,i::Int,::MultiIndexOrderStyle{'C'})=eltype(vs).name.wrapper(indtosub(dims(vs),i,corder),vs)

"""
    searchsortedfirst(vs::AbstractVectorSpace{B},basis::B) where B -> Int
    searchsortedfirst(vs,basis) -> Int

Search the index of a basis in a vector space.
"""
Base.searchsortedfirst(vs::AbstractVectorSpace{B},basis::B) where B=_search_(vs,basis,vs|>typeof|>HasTable)
_search_(vs::AbstractVectorSpace{B},basis::B,::HasTable{true}) where B=_search_(vs,basis,vs|>typeof|>TableSorted,vs|>typeof|>IsMultiIndexable)
_search_(vs::AbstractVectorSpace{B},basis::B,::TableSorted{true},::IsMultiIndexable{false}) where B=searchsortedfirst(vs.table,basis)
_search_(vs::AbstractVectorSpace{B},basis::B,::TableSorted{false},::IsMultiIndexable{false}) where B=findfirst(isequal(basis),vs.table)
_search_(vs::AbstractVectorSpace{B},basis::B,::TableSorted{true},::IsMultiIndexable{true}) where B=searchsortedfirst(vs.table,inds(basis,vs))
_search_(vs::AbstractVectorSpace{B},basis::B,::TableSorted{false},::IsMultiIndexable{true}) where B=findfirst(isequal(inds(basis,vs)),vs.table)
_search_(vs::AbstractVectorSpace{B},basis::B,::HasTable{false}) where B=_search_(vs,basis,vs|>typeof|>IsMultiIndexable)
_search_(vs::AbstractVectorSpace{B},basis::B,::IsMultiIndexable{true}) where B=_search_(vs,basis,vs|>typeof|>MultiIndexOrderStyle)
function _search_(vs::AbstractVectorSpace{B},basis::B,::MultiIndexOrderStyle{'F'}) where B
    basis ∈ vs ? subtoind(dims(vs),inds(basis,vs),forder) : error("searchsortedfirst error: input basis($basis) not found.")
end
function _search_(vs::AbstractVectorSpace{B},basis::B,::MultiIndexOrderStyle{'C'}) where B
    basis ∈ vs ? subtoind(dims(vs),inds(basis,vs),corder) : error("searchsortedfirst error: input basis($basis) not found.")
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
    findfirst(basis::B,vs::AbstractVectorSpace{B}) where B -> Int
    findfirst(bases,vs::AbstractVectorSpace) -> NTuple{length(bases),Int}

Get the index of a basis or the indices of a couple of bases in a vector space.
"""
Base.findfirst(basis::B,vs::AbstractVectorSpace{B}) where B=_find_(basis,vs,vs|>typeof|>TableSorted)
Base.findfirst(bases,vs::AbstractVectorSpace)=NTuple{length(bases),Int}(findfirst(basis,vs)::Int for basis in bases)
_find_(basis::B,vs::AbstractVectorSpace{B},::TableSorted{true}) where B=(index=searchsortedfirst(vs,basis);0<index<=dimension(vs) && vs[index]==basis ? index : missing)
_find_(basis::B,vs::AbstractVectorSpace{B},::TableSorted{false}) where B=searchsortedfirst(vs,basis)

"""
    in(basis::B,vs::AbstractVectorSpace{B}) where B

Judge whether a basis belongs to a vector space.
"""
Base.in(basis::B,vs::AbstractVectorSpace{B}) where B=_in_(basis,vs,vs|>typeof|>HasTable,vs|>typeof|>IsMultiIndexable)
_in_(basis::B,vs::AbstractVectorSpace{B},::HasTable,::IsMultiIndexable) where B=isa(findfirst(basis,vs),Integer)
_in_(basis::B,vs::AbstractVectorSpace{B},::HasTable{false},::IsMultiIndexable{true}) where B=_in_(basis,vs,Val(vs|>typeof|>rank))
@generated function _in_(basis::B,vs::AbstractVectorSpace{B},::Val{N}) where {B,N}
    N==0 && return :(true)
    prepare=:(INDS=inds(basis,vs);DIMS=dims(vs))
    expr=:(0<INDS[1]<=DIMS[1])
    for i=2:N
        expr=Expr(:&&,expr,:(0<INDS[$i]<=DIMS[$i]))
    end
    return Expr(:block,prepare,expr)
end

"""
    DirectVectorSpace{S}(table::NTuple{N,B}) where {S,B,N}
    DirectVectorSpace{S}(table...) where S

Simplest vector space, whose bases are stored in the attribute `:table` as an ntuple.

The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted.
"""
struct DirectVectorSpace{S,B,N} <: AbstractVectorSpace{B}
    table::NTuple{N,B}
    DirectVectorSpace{S}(table::NTuple{N,B}) where {S,B,N}=(@assert(S∈('T','F'),"DirectVectorSpace error: wrong type parameter S($S).");new{S,B,N}(table))
end
DirectVectorSpace{S}(table...) where S=DirectVectorSpace{S}(table)
HasTable(::Type{<:DirectVectorSpace})=HasTable(true)
TableSorted(::Type{<:DirectVectorSpace{'T'}})=TableSorted(true)
TableSorted(::Type{<:DirectVectorSpace{'F'}})=TableSorted(false)

"""
    ⊕(basis1::B,basis2::B) where B
    ⊕(basis::B,vs::DirectVectorSpace{<:Any,B}) where B
    ⊕(vs::DirectVectorSpace{<:Any,B},basis::B) where B
    ⊕(vs1::DirectVectorSpace{<:Any,B},vs2::DirectVectorSpace{<:Any,B}) where B

Get the direct sum between bases or direct vector spaces.
"""
⊕(basis1::B,basis2::B) where B=DirectVectorSpace{'F'}(basis1,basis2)
⊕(basis::B,vs::DirectVectorSpace{<:Any,B}) where B=DirectVectorSpace{'F'}(basis,vs.table...)
⊕(vs::DirectVectorSpace{<:Any,B},basis::B) where B=DirectVectorSpace{'F'}(vs.table...,basis)
⊕(vs1::DirectVectorSpace{<:Any,B},vs2::DirectVectorSpace{<:Any,B}) where B=DirectVectorSpace{'F'}(vs1.table...,vs2.table...)

"""
    AbstractOrderedIndices{N} <: AbstractVectorSpace{NTuple{N,Int}}

A simplest class of multiindexable vector spaces, whose bases are just tuples of integers.

This class of vector spaces must have the following attribute:
`dims::NTuple{N,Int}`: the dimesnions of the Cartesian indices along all axes
"""
abstract type AbstractOrderedIndices{N} <: AbstractVectorSpace{NTuple{N,Int}} end
IsMultiIndexable(::Type{<:AbstractOrderedIndices})=IsMultiIndexable(true)
dims(vs::AbstractOrderedIndices)=vs.dims
inds(basis::NTuple{N,Int},::AbstractOrderedIndices{N}) where N=basis
rank(::Type{<:AbstractOrderedIndices{N}}) where N=N
Tuple(index::NTuple{N,Int},::AbstractOrderedIndices{N}) where N=index

"""
    DirectOrderedIndices{M}(dims::NTuple{N,Int}) where {M,N}

Direct ordered Cartesian indices.

!!! note
    1) It can be C/C++ ordered or Fortran ordered depending on the first type parameter `M`, with `'C'` for the former and `'F'` the latter.
    2) For its bases (Cartesian indices), there is no restriction except that they should be in the proper range defined by its `dims`.
"""
struct DirectOrderedIndices{M,N} <: AbstractOrderedIndices{N}
    dims::NTuple{N,Int}
    DirectOrderedIndices{M}(dims::NTuple{N,Int}) where {M,N}=(@assert(M∈('C','F'),"DirectOrderedIndices error: wrong type parameter M($M).");new{M,N}(dims))
end
DirectOrderedIndices{M}(dims::Int...) where M=DirectOrderedIndices{M}(dims)
MultiIndexOrderStyle(::Type{<:DirectOrderedIndices{M}}) where M=MultiIndexOrderStyle(M)

"""
    TabledOrderedIndices{S}(dims::NTuple{N,Int},table::Vector{NTuple{N,Int}}) where {S,N}
    TabledOrderedIndices{N}(::Type{M},len::Int) where {N,M<:AbstractCombinatorics}

Tabled ordered Cartesian indices.

Compared to [`DirectOrderedIndices`](@ref), the bases of this kind of vector spaces are stored in the attribute `:table`, which must be a vector of tuple of integers. The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted. This type suits the situations when the Cartesian indices are restricted by extra conditions except that propoesed by the attribute `:dims`.
"""
struct TabledOrderedIndices{S,N} <: AbstractOrderedIndices{N}
    dims::NTuple{N,Int}
    table::Vector{NTuple{N,Int}}
    TabledOrderedIndices{S}(dims::NTuple{N,Int},table::Vector{NTuple{N,Int}}) where {S,N}=new{S,N}(dims,table)
end
TabledOrderedIndices{N}(::Type{M},len::Int) where {N,M<:AbstractCombinatorics}=TabledOrderedIndices{'T'}(ntuple(i->len,N),M{N}(1:len)|>collect)
HasTable(::Type{<:TabledOrderedIndices})=HasTable(true)
TableSorted(::Type{<:TabledOrderedIndices{'T'}})=TableSorted(true)
TableSorted(::Type{<:TabledOrderedIndices{'F'}})=TableSorted(false)

end # module
