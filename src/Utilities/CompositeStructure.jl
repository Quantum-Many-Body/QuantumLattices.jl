module CompositeStructure

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
function Base.:(==)(ct1::CompositeNTuple,ct2::CompositeNTuple)
    fieldcount(ct1|>typeof)==fieldcount(ct2|>typeof) ? all(getfield(ct1,i)==getfield(ct2,i) for i=1:fieldcount(ct1|>typeof)) : false
end
function Base.isequal(ct1::CompositeNTuple,ct2::CompositeNTuple)
    fieldcount(ct1|>typeof)==fieldcount(ct2|>typeof) ? all(isequal(getfield(ct1,i),getfield(ct2,i)) for i=1:fieldcount(ct1|>typeof)) : false
end
Base.getindex(ct::CompositeNTuple,i::Union{<:Integer,CartesianIndex})=ct.contents[i]
Base.getindex(ct::CompositeNTuple,inds)=typeof(ct).name.wrapper((name==:contents ? ct.contents[inds] : getfield(ct,name) for name in ct|>typeof|>fieldnames)...)
Base.iterate(ct::CompositeNTuple)=iterate(ct.contents)
Base.iterate(ct::CompositeNTuple,state)=iterate(ct.contents,state)
Base.keys(ct::CompositeNTuple)=keys(ct.contents)
Base.values(ct::CompositeNTuple)=values(ct.contents)
Base.pairs(ct::CompositeNTuple)=pairs(ct.contents)
Base.convert(::Type{Tuple},ct::CompositeNTuple)=ct.contents

"""
    CompositeVector{T}

A composite vector is a vector that is implemented by including an ordinary `Vector` as one of its attributes with the name `:contents`.
"""
abstract type CompositeVector{T} <:AbstractVector{T} end
Base.size(cv::CompositeVector)=size(cv.contents)
Base.size(cv::CompositeVector,i)=size(cv.contents,i)
Base.length(cv::CompositeVector)=length(cv.contents)
function Base.:(==)(cv1::CompositeVector,cv2::CompositeVector)
    fieldcount(cv1|>typeof)==fieldcount(cv2|>typeof) ? all(getfield(cv1,i)==getfield(cv2,i) for i=1:fieldcount(cv1|>typeof)) : false
end
function Base.isequal(cv1::CompositeVector,cv2::CompositeVector)
    fieldcount(cv1|>typeof)==fieldcount(cv2|>typeof) ? all(isequal(getfield(cv1,i),getfield(cv2,i)) for i=1:fieldcount(cv1|>typeof)) : false
end
Base.getindex(cv::CompositeVector,i::Union{<:Integer,CartesianIndex})=cv.contents[i]
Base.getindex(cv::CompositeVector,inds)=typeof(cv)((name==:contents ? cv.contents[inds] : getfield(cv,name) for name in cv|>typeof|>fieldnames)...)
Base.setindex!(cv::CompositeVector,value,inds)=(cv.contents[inds]=value)
Base.push!(cv::CompositeVector,values...)=(push!(cv.contents,values...);cv)
Base.pushfirst!(cv::CompositeVector,values...)=(pushfirst!(cv.contents,values...);cv)
Base.insert!(cv::CompositeVector,index::Integer,value)=(insert!(cv.contents,index,value);cv)
Base.append!(cv::CompositeVector,values)=(append!(cv.contents,values);cv)
Base.prepend!(cv::CompositeVector,values)=(prepend!(cv.contents,values);cv)
Base.splice!(cv::CompositeVector,index::Integer,replacement=Base._default_splice)=splice!(cv.contents,index,replacement)
function Base.splice!(cv::CompositeVector,range::UnitRange{<:Integer},replacement=Base._default_splice)
    rmvs=splice!(cv.contents,range,replacement)
    typeof(cv)((name==:contents ? rmvs : getfield(cv,name) for name in cv|>typeof|>fieldnames)...)
end
Base.deleteat!(cv::CompositeVector,indices)=(deleteat!(cv.contents,indices);cv)
Base.pop!(cv::CompositeVector)=pop!(cv.contents)
Base.popfirst!(cv::CompositeVector)=popfirst!(cv.contents)
Base.empty!(cv::CompositeVector)=(empty!(cv.contents);cv)
Base.empty(cv::CompositeVector)=typeof(cv)((name==:contents ? empty(cv.contents) : getfield(cv,name) for name in cv|>typeof|>fieldnames)...)
Base.reverse(cv::CompositeVector)=typeof(cv)((name==:contents ? reverse(cv.contents) : getfield(cv,name) for name in cv|>typeof|>fieldnames)...)
function Base.similar(cv::CompositeVector,dtype::Type{T}=eltype(cv),dims::NTuple{N,Int}=size(cv)) where {T,N}
    typeof(cv).name.wrapper((name==:contents ? similar(cv.contents,dtype,dims) : getfield(cv,name) for name in cv|>typeof|>fieldnames)...)
end
Base.iterate(cv::CompositeVector,state=1)=iterate(cv.contents,state)
Base.keys(cv::CompositeVector)=keys(cv.contents)
Base.values(cv::CompositeVector)=values(cv.contents)
Base.pairs(cv::CompositeVector)=pairs(cv.contents)
Base.convert(::Type{Vector},cv::CompositeVector)=cv.contents

"""
    CompositeDict{K,V}

A composite dict is a dict that is implemented by including an ordinary `Dict` as one of its attributes with the name `:contents`.
"""
abstract type CompositeDict{K,V} <: AbstractDict{K,V} end
Base.isempty(cd::CompositeDict)=isempty(cd.contents)
Base.length(cd::CompositeDict)=length(cd.contents)
Base.haskey(cd::CompositeDict,key)=haskey(cd.contents,key)
Base.in(p::Pair,cd::CompositeDict,valcmp=(==))=in(p,cd.contents,valcmp)
Base.hash(cd::CompositeDict,h::UInt)=hash(Tuple(getfield(cd,name) for name in cd|>typeof|>fieldnames),h)
function Base.:(==)(cd1::CompositeDict,cd2::CompositeDict)
    fieldcount(cd1|>typeof)==fieldcount(cd2|>typeof) ? all(getfield(cd1,i)==getfield(cd2,i) for i=1:fieldcount(cd1|>typeof)) : false
end
function Base.isequal(cd1::CompositeDict,cd2::CompositeDict)
    fieldcount(cd1|>typeof)==fieldcount(cd2|>typeof) ? all(isequal(getfield(cd1,i),getfield(cd2,i)) for i=1:fieldcount(cd1|>typeof)) : false
end
Base.get(cd::CompositeDict,key,default)=get(cd.contents,key,default)
Base.get(f::Base.Callable,cd::CompositeDict,key)=get(f,cd.contents,key)
Base.getkey(cd::CompositeDict,key,default)=getkey(cd.contents,key,default)
Base.getindex(cd::CompositeDict{K,V},index::K) where {K,V}=cd.contents[index]
Base.push!(cd::CompositeDict,ps::Pair...)=(push!(cd.contents,ps...);cd)
Base.get!(cd::CompositeDict,key,default)=get!(cd.contents,key,default)
Base.get!(default,cd::CompositeDict,key)=get!(default,cd.contents,key)
Base.setindex!(cd::CompositeDict{K,V},value::V,index::K) where {K,V}=(cd.contents[index]=value)
Base.pop!(cd::CompositeDict)=pop!(cd.contents)
Base.pop!(cd::CompositeDict,key)=pop!(cd.contents,key)
Base.pop!(cd::CompositeDict,key,default)=pop!(cd.contents,key,default)
Base.delete!(cd::CompositeDict,key)=(delete!(cd.contents,key);cd)
Base.empty!(cd::CompositeDict)=(empty!(cd.contents);cd)
Base.merge(cd::CD,others::CD...) where CD<:CompositeDict=merge!(empty(cd),cd,others...)
Base.merge(combine::Function,cd::CD,others::CD...) where CD<:CompositeDict=merge!(combine,empty(cd),cd,others...)
Base.empty(cd::CompositeDict)=typeof(cd)((name==:contents ? Dict{cd|>keytype,cd|>valtype}() : getfield(cd,name) for name in fieldnames(typeof(cd)))...)
Base.iterate(cd::CompositeDict)=iterate(cd.contents)
Base.iterate(cd::CompositeDict,state)=iterate(cd.contents,state)
Base.keys(cd::CompositeDict)=keys(cd.contents)
Base.values(cd::CompositeDict)=values(cd.contents)
Base.pairs(cd::CompositeDict)=pairs(cd.contents)
Base.convert(::Type{Dict},cd::CompositeDict)=cd.contents

end #module
