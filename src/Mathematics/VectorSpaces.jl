module VectorSpaces

using ..Combinatorics: AbstractCombinatorics
using ...Prerequisites.TypeTraits: efficientoperations,forder,corder,subtoind,indtosub

import ...Interfaces: dimension,⊕,rank,dims,inds,degree

export VectorSpace
export HasTable,TableSorted,IsMultiIndexable,MultiIndexOrderStyle
export SimpleVectorSpace
export OrderedIndices,SimpleIndices,TabledIndices
export GradedTables,GradedVectorSpace

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
_getindex_(vs::VectorSpace,i::Int,::HasTable{true},::IsMultiIndexable{true})=eltype(vs).name.wrapper(vs.table[i],vs)
_getindex_(vs::VectorSpace,i::Int,::HasTable{false},::IsMultiIndexable{true})=_getindex_(vs,i,vs|>typeof|>MultiIndexOrderStyle)
_getindex_(vs::VectorSpace,i::Int,::MultiIndexOrderStyle{'F'})=eltype(vs).name.wrapper(indtosub(dims(vs),i,forder),vs)
_getindex_(vs::VectorSpace,i::Int,::MultiIndexOrderStyle{'C'})=eltype(vs).name.wrapper(indtosub(dims(vs),i,corder),vs)

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
    GradedTables(vs::Tuple,ks::Tuple)
    GradedTables(::Type{M},n::Int,gs::Val{grades}) where {M<:AbstractCombinatorics,grades}

The tables of a graded vector space.

Alias for `Base.Iterators.Pairs{G,V,KS<:Tuple,VS<:Tuple}`.
"""
const GradedTables{G,V,KS<:Tuple,VS<:Tuple}=Base.Iterators.Pairs{G,V,KS,VS}
GradedTables(vs::Tuple,ks::Tuple)=Base.Iterators.Pairs(vs,ks)
function GradedTables(::Type{M},n::Int,gs::Val{grades}) where {M<:AbstractCombinatorics,grades}
    @assert isa(grades,Tuple{Vararg{Int}}) "GradedTables error: grades should be tuple of integers."
    return GradedTables(comdata(M,n,gs),grades)
end
@generated function comdata(::Type{M},n::Int,::Val{grades}) where {M<:AbstractCombinatorics,grades}
    return Expr(:tuple,[:(TabledIndices{$g}(M,n)) for g in grades]...)
end

"""
    keytype(tables::GradedTables,g::Int)
    keytype(::Type{T},g::Int) where {T<:GradedTables}

Get the gth keytype of a graded tables.
"""
Base.keytype(tables::GradedTables,g::Int)=keytype(typeof(tables),g)
Base.keytype(::Type{T},g::Int) where {T<:GradedTables}=fieldtype(fieldtype(T,:itr),g)

"""
    valtype(tables::GradedTables,g::Int)
    valtype(::Type{T},g::Int) where {T<:GradedTables}

Get the gth valtype of a graded tables.
"""
Base.valtype(tables::GradedTables,g::Int)=valtype(typeof(tables),g)
Base.valtype(::Type{T},g::Int) where {T<:GradedTables}=fieldtype(fieldtype(T,:data),g)

"""
    rank(tables::GradedTables) -> Int
    rank(::Type{T}) where {T<:GradedTables} -> Int

Get the total number of keys or values of a graded tables.
"""
rank(tables::GradedTables)=tables|>typeof|>rank
rank(::Type{T}) where {T<:GradedTables}=fieldcount(fieldtype(T,:data))

"""
    GradedVectorSpace{G,B,V<:VectorSpace,T<:GradedTables{G,V}} <: VectorSpace{Tuple{G,B}}

Abstract type of graded vector spaces.

A graded vector space is a vector space that has the extra structure of a grading, which is a decomposition of the vector space into a direct sum of vector subspaces.

Concrete subtypes must have the following attribute:
* `:tables::GradedTables{G,V} where {G,V<:VectorSpace}`: the tables of the subspaces.
"""
abstract type GradedVectorSpace{G,B,V<:VectorSpace,T<:GradedTables{G,V}} <: VectorSpace{Tuple{G,B}} end

"""
    keys(vs::GradedVectorSpace) -> Tuple
    values(vs::GradedVectorSpace) -> Tuple
    pairs(vs::GradedVectorSpace) -> GradedTables

Iterate over the keys, values or pairs of a graded vector space.
"""
Base.keys(vs::GradedVectorSpace)=keys(vs.tables)
Base.values(vs::GradedVectorSpace)=values(vs.tables)
Base.pairs(vs::GradedVectorSpace)=pairs(vs.tables)

"""
    keytype(vs::GradedVectorSpace,g::Int)
    keytype(::Type{V},g::Int) where {V<:GradedVectorSpace}

Get the gth keytype of a graded vector space.
"""
Base.keytype(vs::GradedVectorSpace,g::Int)=keytype(typeof(vs),g)
Base.keytype(::Type{V},g::Int) where {V<:GradedVectorSpace}=keytype(fieldtype(V,:tables),g)

"""
    valtype(vs::GradedVectorSpace,g::Int)
    valtype(::Type{V},g::Int) where {V<:GradedVectorSpace}

Get the gth valtype of a graded vector space.
"""
Base.valtype(vs::GradedVectorSpace,g::Int)=valtype(typeof(vs),g)
Base.valtype(::Type{V},g::Int) where {V<:GradedVectorSpace}=valtype(fieldtype(V,:tables),g)

"""
    eltype(vs::GradedVectorSpace,g::Int)
    eltype(::Type{V},g::Int) where {V<:GradedVectorSpace}

Get the gth eltype of a graded vector space.
"""
Base.eltype(vs::GradedVectorSpace,g::Int)=eltype(typeof(vs),g)
Base.eltype(::Type{V},g::Int) where {V<:GradedVectorSpace}=Tuple{keytype(V,g),eltype(valtype(V,g))}

"""
    rank(vs::GradedVectorSpace) -> Int
    rank(::Type{V}) where {V<:GradedVectorSpace} -> Int

Get the rank, i.e. the total number of vector subspaces.
"""
rank(vs::GradedVectorSpace)=vs|>typeof|>rank
rank(::Type{V}) where {V<:GradedVectorSpace}=rank(fieldtype(V,:tables))

"""
    degree(g::G,vs::GradedVectorSpace{G}) where G -> Int

Get the degree of a vector subspace whose grade are represented by g.
"""
degree(g::G,vs::GradedVectorSpace{G}) where G=findfirst(isequal(g),keys(vs))

"""
    dimension(vs::GradedVectorSpace) -> Int
    dimension(vs::GradedVectorSpace{G},g::G) where G -> Int
    dimension(vs::GradedVectorSpace{G},gs::NTuple{N,G}) where {G,N} -> Int

Get the dimension of the whole graded vector space or some vector subspaces.
"""
dimension(vs::GradedVectorSpace)=sum(NTuple{rank(typeof(vs)),Int}(dimension(v) for v in values(vs)))::Int
dimension(vs::GradedVectorSpace{G},g::G) where G=(index=degree(g,vs);isa(index,Int) ? dimension(values(vs)[index]) : 0)::Int
function dimension(vs::GradedVectorSpace{G},gs::NTuple{N,G}) where {G,N}
    result=0
    for g in gs
        result=result+dimension(vs,g)::Int
    end
    result
end

"""
    iterate(vs::GradedVectorSpace,state=(1,1))

Iterate over the whole graded vector space.
"""
function Base.iterate(vs::GradedVectorSpace,state=(1,1))
    tables,degree,index=values(vs),state[1],state[2]
    while index>length(tables[degree])
        degree=degree+1
        degree>length(tables) && return nothing
        index=1
    end
    return (keys(vs)[degree],tables[degree][index]),(degree,index+1)
end

"""
    getindex(vs::GradedVectorSpace{B},pair::Tuple{B,Int}) where B
    getindex(vs::GradedVectorSpace,i::Int)

Get the basis of a graded vector space by a grade-index pair or by an index.
"""
Base.getindex(vs::GradedVectorSpace{B},pair::Tuple{B,Int}) where B=values(vs)[degree(pair[1],vs)][pair[2]]
function Base.getindex(vs::GradedVectorSpace,i::Int)
    degree,dimensions=1,NTuple{rank(typeof(vs)),Int}(dimension(v) for v in values(vs))
    while i>dimensions[degree]
        i=i-dimensions[degree]
        degree=degree+1
        degree>length(dimensions) && error("getindex error: input index($i) out of range.")
    end
    return (keys(vs)[degree],values(vs)[degree][i])
end

"""
    searchsortedfirst(vs::GradedVectorSpace{G,B},pair::Tuple{G,B}) where {G,B} -> Int

Find the index of a grade-basis pair in a graded vector space.
"""
function Base.searchsortedfirst(vs::GradedVectorSpace{G,B},pair::Tuple{G,B}) where {G,B}
    dim,tables,index,basis=0,values(vs),degree(pair[1],vs),pair[2]
    for i=1:(index-1)
        dim=dim+dimension(tables[i])
    end
    searchsortedfirst(tables[index],basis)+dim
end

end # module
