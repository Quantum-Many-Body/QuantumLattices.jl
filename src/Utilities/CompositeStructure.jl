module CompositeStructure

using ..Utilities: comparison

export CompositeNTuple,CompositeVector,CompositeDict

"""
    CompositeNTuple{N,T}

A composite ntuple is a ntuple that is implemented by including an ordinary `NTuple` as one of its attributes with the name `:contents`.
"""
abstract type CompositeNTuple{N,T} end
Base.length(::CompositeNTuple{N,T}) where {N,T}=N
Base.length(::Type{<:CompositeNTuple{N,T}}) where {N,T}=N
Base.eltype(::CompositeNTuple{N,T}) where {N,T}=T
Base.eltype(::Type{<:CompositeNTuple{N,T}}) where {N,T}=T
Base.:(==)(ct1::CompositeNTuple,ct2::CompositeNTuple) = ==(comparison,ct1,ct2)
Base.isequal(ct1::CompositeNTuple,ct2::CompositeNTuple)=isequal(comparison,ct1,ct2)
Base.getindex(ct::CompositeNTuple,i::Union{<:Integer,CartesianIndex})=getfield(ct,:contents)[i]
@generated function Base.getindex(ct::CompositeNTuple,inds)
    exprs=[name==:contents ? :(getfield(ct,:contents)[inds]) : :(getfield(ct,$i)) for (i,name) in enumerate(ct|>fieldnames)]
    return :(typeof(ct).name.wrapper($(exprs...)))
end
Base.iterate(ct::CompositeNTuple)=iterate(getfield(ct,:contents))
Base.iterate(ct::CompositeNTuple,state)=iterate(getfield(ct,:contents),state)
Base.keys(ct::CompositeNTuple)=keys(getfield(ct,:contents))
Base.values(ct::CompositeNTuple)=values(getfield(ct,:contents))
Base.pairs(ct::CompositeNTuple)=pairs(getfield(ct,:contents))
Base.convert(::Type{Tuple},ct::CompositeNTuple)=getfield(ct,:contents)

"""
    CompositeVector{T}

A composite vector is a vector that is implemented by including an ordinary `Vector` as one of its attributes with the name `:contents`.
"""
abstract type CompositeVector{T} <:AbstractVector{T} end
Base.size(cv::CompositeVector)=size(getfield(cv,:contents))
Base.size(cv::CompositeVector,i)=size(getfield(cv,:contents),i)
Base.length(cv::CompositeVector)=length(getfield(cv,:contents))
Base.:(==)(cv1::CompositeVector,cv2::CompositeVector) = ==(comparison,cv1,cv2)
Base.isequal(cv1::CompositeVector,cv2::CompositeVector)=isequal(comparison,cv1,cv2)
Base.getindex(cv::CompositeVector,i::Union{<:Integer,CartesianIndex})=getfield(cv,:contents)[i]
@generated function Base.getindex(cv::CompositeVector,inds)
    exprs=[name==:contents ? :(getfield(cv,:contents)[inds]) : :(getfield(cv,$i)) for (i,name) in enumerate(cv|>fieldnames)]
    return :(typeof(cv)($(exprs...)))
end
Base.setindex!(cv::CompositeVector,value,inds)=(getfield(cv,:contents)[inds]=value)
Base.push!(cv::CompositeVector,values...)=(push!(getfield(cv,:contents),values...);cv)
Base.pushfirst!(cv::CompositeVector,values...)=(pushfirst!(getfield(cv,:contents),values...);cv)
Base.insert!(cv::CompositeVector,index::Integer,value)=(insert!(getfield(cv,:contents),index,value);cv)
Base.append!(cv::CompositeVector,values)=(append!(getfield(cv,:contents),values);cv)
Base.prepend!(cv::CompositeVector,values)=(prepend!(getfield(cv,:contents),values);cv)
Base.splice!(cv::CompositeVector,index::Integer,replacement=Base._default_splice)=splice!(getfield(cv,:contents),index,replacement)
@generated function Base.splice!(cv::CompositeVector,range::UnitRange{<:Integer},replacement=Base._default_splice)
    exprs=[name==:contents ? :(splice!(getfield(cv,:contents),range,replacement)) : :(getfield(cv,$i)) for (i,name) in enumerate(cv|>fieldnames)]
    return :(typeof(cv)($(exprs...)))
end
Base.deleteat!(cv::CompositeVector,indices)=(deleteat!(getfield(cv,:contents),indices);cv)
Base.pop!(cv::CompositeVector)=pop!(getfield(cv,:contents))
Base.popfirst!(cv::CompositeVector)=popfirst!(getfield(cv,:contents))
Base.empty!(cv::CompositeVector)=(empty!(getfield(cv,:contents));cv)
@generated function Base.empty(cv::CompositeVector)
    exprs=[name==:contents ? :(empty(getfield(cv,:contents))) : :(getfield(cv,$i)) for (i,name) in enumerate(cv|>fieldnames)]
    return :(typeof(cv)($(exprs...)))
end
@generated function Base.reverse(cv::CompositeVector)
    exprs=[name==:contents ? :(reverse(getfield(cv,:contents))) : :(getfield(cv,$i)) for (i,name) in enumerate(cv|>fieldnames)]
    return :(typeof(cv)($(exprs...)))
end
@generated function Base.similar(cv::CompositeVector,dtype::Type{T}=eltype(cv),dims::NTuple{N,Int}=size(cv)) where {T,N}
    exprs=[name==:contents ? :(similar(getfield(cv,:contents),dtype,dims)) : :(getfield(cv,$i)) for (i,name) in enumerate(cv|>fieldnames)]
    return :(typeof(cv).name.wrapper($(exprs...)))
end
Base.iterate(cv::CompositeVector,state=1)=iterate(getfield(cv,:contents),state)
Base.keys(cv::CompositeVector)=keys(getfield(cv,:contents))
Base.values(cv::CompositeVector)=values(getfield(cv,:contents))
Base.pairs(cv::CompositeVector)=pairs(getfield(cv,:contents))
Base.convert(::Type{Vector},cv::CompositeVector)=getfield(cv,:contents)

"""
    CompositeDict{K,V}

A composite dict is a dict that is implemented by including an ordinary `Dict` as one of its attributes with the name `:contents`.
"""
abstract type CompositeDict{K,V} <: AbstractDict{K,V} end
Base.isempty(cd::CompositeDict)=isempty(getfield(cd,:contents))
Base.length(cd::CompositeDict)=length(getfield(cd,:contents))
Base.haskey(cd::CompositeDict,key)=haskey(getfield(cd,:contents),key)
Base.in(p::Pair,cd::CompositeDict,valcmp=(==))=in(p,getfield(cd,:contents),valcmp)
Base.hash(cd::CompositeDict,h::UInt)=hash(Tuple(getfield(cd,name) for name in cd|>typeof|>fieldnames),h)
Base.:(==)(cd1::CompositeDict,cd2::CompositeDict) = ==(comparison,cd1,cd2)
Base.isequal(cd1::CompositeDict,cd2::CompositeDict)=isequal(comparison,cd1,cd2)
Base.get(cd::CompositeDict,key,default)=get(getfield(cd,:contents),key,default)
Base.get(f::Base.Callable,cd::CompositeDict,key)=get(f,getfield(cd,:contents),key)
Base.getkey(cd::CompositeDict,key,default)=getkey(getfield(cd,:contents),key,default)
Base.getindex(cd::CompositeDict{K,V},index::K) where {K,V}=getfield(cd,:contents)[index]
Base.push!(cd::CompositeDict,ps::Pair...)=(push!(getfield(cd,:contents),ps...);cd)
Base.get!(cd::CompositeDict,key,default)=get!(getfield(cd,:contents),key,default)
Base.get!(default,cd::CompositeDict,key)=get!(default,getfield(cd,:contents),key)
Base.setindex!(cd::CompositeDict{K,V},value::V,index::K) where {K,V}=(getfield(cd,:contents)[index]=value)
Base.pop!(cd::CompositeDict)=pop!(getfield(cd,:contents))
Base.pop!(cd::CompositeDict,key)=pop!(getfield(cd,:contents),key)
Base.pop!(cd::CompositeDict,key,default)=pop!(getfield(cd,:contents),key,default)
Base.delete!(cd::CompositeDict,key)=(delete!(getfield(cd,:contents),key);cd)
Base.empty!(cd::CompositeDict)=(empty!(getfield(cd,:contents));cd)
Base.merge(cd::CD,others::CD...) where CD<:CompositeDict=merge!(empty(cd),cd,others...)
Base.merge(combine::Function,cd::CD,others::CD...) where CD<:CompositeDict=merge!(combine,empty(cd),cd,others...)
@generated function Base.empty(cd::CompositeDict)
    exprs=[name==:contents ? :(Dict{cd|>keytype,cd|>valtype}()) : :(getfield(cd,$i)) for (i,name) in enumerate(cd|>fieldnames)]
    return :(typeof(cd)($(exprs...)))
end
Base.iterate(cd::CompositeDict)=iterate(getfield(cd,:contents))
Base.iterate(cd::CompositeDict,state)=iterate(getfield(cd,:contents),state)
Base.keys(cd::CompositeDict)=keys(getfield(cd,:contents))
Base.values(cd::CompositeDict)=values(getfield(cd,:contents))
Base.pairs(cd::CompositeDict)=pairs(getfield(cd,:contents))
Base.convert(::Type{Dict},cd::CompositeDict)=getfield(cd,:contents)

end #module
