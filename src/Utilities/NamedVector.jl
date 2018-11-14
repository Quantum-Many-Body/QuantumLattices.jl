module NamedVector

import Printf: @printf

export AbstractNamedVector
export @namedvector

"""
    AbstractNamedVector{T}

Abstract type for all concrete named vectors.
"""
abstract type AbstractNamedVector{T} end

"""
    getindex(nv::AbstractNamedVector,index::Int)

Get the value by the `[]` syntax.
"""
Base.getindex(nv::AbstractNamedVector,index::Int)=getfield(nv,index)

"""
    setindex!(nv::AbstractNamedVector,value,index::Int)

Set the value by the `[]` syntax if mutable.
"""
Base.setindex!(nv::AbstractNamedVector,value,index::Int)=setfield!(nv,index,value)

"""
    ==(nv1::AbstractNamedVector,nv2::AbstractNamedVector) -> Bool

Overloaded `==` operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.
!!! note
    It is not necessary for two named vectors to be of the same concrete type to be equal to each other.
"""
Base.:(==)(nv1::AbstractNamedVector,nv2::AbstractNamedVector)=nv1|>keys==nv2|>keys && nv1|>values==nv2|>values

"""
    <(nv1:NV,nv2:NV) where NV<:AbstractNamedVector -> Bool

Overloaded `<` operator.
"""
Base.:<(nv1::NV,nv2::NV) where NV<:AbstractNamedVector=isless(nv1,nv2)

"""
    isless(nv1::NV,nv2::NV) where NV<:AbstractNamedVector -> Bool

Overloaded `isless` function.
"""
@generated function Base.isless(nv1::NV,nv2::NV) where NV<:AbstractNamedVector
    N=NV|>fieldnames|>length
    expr=Expr(:if,:(getfield(nv1,$N)<getfield(nv2,$N)),true,false)
    for i in range(N-1,stop=1,step=-1)
        expr=Expr(:if,:(getfield(nv1,$i)>getfield(nv2,$i)),false,expr)
        expr=Expr(:if,:(getfield(nv1,$i)<getfield(nv2,$i)),true,expr)
    end
    return expr
end

"""
    show(io::IO,nv::AbstractNamedVector)

Show a concrete `AbstractNamedVector`.
"""
Base.show(io::IO,nv::AbstractNamedVector)=@printf io "%s(%s)" nv|>typeof|>nameof join(nv|>values,',')

"""
    hash(nv::AbstractNamedVector,h::UInt)

Hash a concrete `AbstractNamedVector`.
"""
Base.hash(nv::AbstractNamedVector,h::UInt)=hash(nv|>values,h)

"""
    convert(::Type{Tuple},nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}
    convert(::Type{NTuple},nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}
    convert(::Type{NTuple{N,T}},nv::AbstractNamedVector{T}) where {N,T} -> NTuple{nv|>length,nv|>eltype}

Convert a named vector to tuple.
"""
Base.convert(::Type{Tuple},nv::AbstractNamedVector)=NTuple{nv|>length,nv|>eltype}(getfield(nv,i) for i=1:length(nv))
Base.convert(::Type{NTuple},nv::AbstractNamedVector)=convert(Tuple,nv)
Base.convert(::Type{NTuple{N,T}},nv::AbstractNamedVector{T}) where {N,T}=convert(Tuple,nv)

"""
    length(::Type{NV}) where NV<:AbstractNamedVector -> Int
    length(nv::AbstractNamedVector) -> Int

Get the length of a concrete `AbstractNamedVector`.
"""
Base.length(::Type{NV}) where NV<:AbstractNamedVector=NV|>fieldnames|>length
Base.length(nv::AbstractNamedVector)=nv|>typeof|>length

"""
    eltype(::Type{NV}) where NV<:AbstractNamedVector{T} where T
    eltype(nv::AbstractNamedVector)

Get the type parameter of a concrete `AbstractNamedVector`.
"""
Base.eltype(::Type{<:AbstractNamedVector{T}}) where T=T
Base.eltype(nv::AbstractNamedVector)=nv|>typeof|>eltype

"""
    zero(::Type{NV}) where NV<:AbstractNamedVector
    zero(nv::AbstractNamedVector)

Get a concrete `AbstractNamedVector` with all values being zero.
"""
@generated Base.zero(::Type{NV}) where NV<:AbstractNamedVector=(zeros=(zero(NV|>eltype) for i=1:length(NV|>fieldnames));:(NV($(zeros...))))
Base.zero(nv::AbstractNamedVector)=nv|>typeof|>zero

"""
    iterate(nv::AbstractNamedVector,state=1)
    iterate(rv::Iterators.Reverse{<:AbstractNamedVector},state=length(rv.itr))

Iterate or reversely iterate over the values of a concrete `AbstractNamedVector`.
"""
Base.iterate(nv::AbstractNamedVector,state=1)=state>length(nv) ? nothing : (getfield(nv,state),state+1)
Base.iterate(rv::Iterators.Reverse{<:AbstractNamedVector},state=length(rv.itr))=state<1 ? nothing : (getfield(rv.itr,state),state-1)

"""
    keys(nv::AbstractNamedVector) -> NTuple(nv|>length,Symbol)

Iterate over the names.
"""
Base.keys(nv::AbstractNamedVector)=nv|>typeof|>fieldnames

"""
    values(nv::AbstractNamedVector) -> NTuple{nv|>length,nv|>eltype}

Iterate over the values.
"""
Base.values(nv::AbstractNamedVector)=NTuple{nv|>length,nv|>typeof|>eltype}(getfield(nv,i) for i=1:length(nv))

"""
    pairs(nv::AbstractNamedVector)

Iterate over the `name=>value` pairs.
"""
Base.pairs(nv::AbstractNamedVector)=Base.Generator(=>,keys(nv),values(nv))

"""
    replace(nv::AbstractNamedVector;kwargs...) -> typeof(nv)

Return a copy of a concrete `AbstractNamedVector` with some of the field values replaced by the keyword arguments.
"""
Base.replace(nv::AbstractNamedVector;kwargs...)=(nv|>typeof)((get(kwargs,key,getfield(nv,key)) for key in nv|>keys)...)

"""
    map(f,nvs::NV...) where NV<:AbstractNamedVector -> NV

Apply function `f` elementwise on the input named vectors.
"""
@generated function Base.map(f,nvs::NV...) where NV<:AbstractNamedVector
    exprs=Vector{Expr}(undef,NV|>fieldnames|>length)
    for i=1:length(exprs)
        tmp=[:(nvs[$j][$i]) for j=1:length(nvs)]
        exprs[i]=:(f($(tmp...)))
    end
    return :(($NV)($(exprs...)))
end

"""
    (::Type{NV})(values::NTuple{N,T}) where {NV<:AbstractNamedVector,N,T}

Construct a concrete named vector by a tuple.
"""
(::Type{NV})(values::NTuple{N,T}) where {NV<:AbstractNamedVector,N,T}=NV(values...)

"""
    @namedvector mutableornot::Bool typename fieldnames dtype::Union{Expr,Symbol}=:nothing supertypename=:AbstractNamedVector

Construct a mutable or immutable concrete named vector with the type name being `typename` and the fieldnames specified by `fieldnames`, and optionally, the type parameters specified by `dtype` and the supertype specified by `supertypename`.
"""
macro namedvector(mutableornot::Bool,typename,fieldnames,dtype::Union{Expr,Symbol}=:nothing,supertypename=:AbstractNamedVector)
    typename,supertypename=Symbol(typename),Symbol(supertypename)
    fieldnames=tuple(eval(fieldnames)...)
    @assert all(isa(name,Symbol) for name in fieldnames) "namedvector error: every field name should be a `Symbol`."
    isa(dtype,Expr) && (@assert (dtype.head==:(<:) && dtype.args|>length==1) "namedvector error: not supported `dtype`.")
    dname,dscope=isa(dtype,Expr) ? (:T,dtype.args[1]) : (dtype==:nothing ? (:T,:Any) : (dtype,:concrete))
    if dscope==:concrete
        new=:($(esc(typename)))
        super=supertypename==:AbstractNamedVector ? :(AbstractNamedVector{$(esc(dname))}) : :($(esc(supertypename)){$(esc(dname))})
        body=(:($field::$(esc(dname))) for field in fieldnames)
    else
        new=:($(esc(typename)){$dname<:$(esc(dscope))})
        super=supertypename==:AbstractNamedVector ? :(AbstractNamedVector{$dname}) : :($(esc(supertypename)){$dname})
        body=(:($field::$dname) for field in fieldnames)
    end
    structdef=Expr(:struct,mutableornot,Expr(:<:,new,super),Expr(:block,body...))
    functions=:(Base.fieldnames(::Type{<:$(esc(typename))})=$fieldnames)
    return Expr(:block,:(Base.@__doc__($structdef)),functions)
end

end #module
