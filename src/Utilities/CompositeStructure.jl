module CompositeStructure

export CompositeVector,CompositeDict

"""
    CompositeVector{T}

A composite vector is a vector that is implemented by including an ordinary `Vector` as one of its attributes with the name `:contents`.
"""
abstract type CompositeVector{T} <:AbstractVector{T} end
Base.size(fv::CompositeVector)=size(fv.contents)
Base.size(fv::CompositeVector,i)=size(fv.contents,i)
Base.length(fv::CompositeVector)=length(fv.contents)
Base.:(==)(fv1::FV,fv2::FV) where FV<:CompositeVector=all(getfield(fv1,name)==getfield(fv2,name) for name in FV|>fieldnames)
Base.isequal(fv1::FV,fv2::FV) where FV<:CompositeVector=all(isequal(getfield(fv1,name),getfield(fv2,name)) for name in FV|>fieldnames)
Base.getindex(fv::CompositeVector,i::Union{<:Integer,CartesianIndex})=fv.contents[i]
Base.getindex(fv::CompositeVector,inds)=typeof(fv)((name==:contents ? fv.contents[inds] : getfield(fv,name) for name in fv|>typeof|>fieldnames)...)
Base.setindex!(fv::CompositeVector,value,inds)=(fv.contents[inds]=value)
Base.push!(fv::CompositeVector,values...)=(push!(fv.contents,values...);fv)
Base.pushfirst!(fv::CompositeVector,values...)=(pushfirst!(fv.contents,values...);fv)
Base.insert!(fv::CompositeVector,index::Integer,value)=(insert!(fv.contents,index,value);fv)
Base.append!(fv::CompositeVector,values)=(append!(fv.contents,values);fv)
Base.prepend!(fv::CompositeVector,values)=(prepend!(fv.contents,values);fv)
Base.splice!(fv::CompositeVector,index::Integer,replacement=Base._default_splice)=splice!(fv.contents,index,replacement)
function Base.splice!(fv::CompositeVector,range::UnitRange{<:Integer},replacement=Base._default_splice)
    rmvs=splice!(fv.contents,range,replacement)
    typeof(fv)((name==:contents ? rmvs : getfield(fv,name) for name in fv|>typeof|>fieldnames)...)
end
Base.deleteat!(fv::CompositeVector,indices)=(deleteat!(fv.contents,indices);fv)
Base.pop!(fv::CompositeVector)=pop!(fv.contents)
Base.popfirst!(fv::CompositeVector)=popfirst!(fv.contents)
Base.empty!(fv::CompositeVector)=(empty!(fv.contents);fv)
Base.empty(fv::CompositeVector)=typeof(fv)((name==:contents ? empty(fv.contents) : getfield(fv,name) for name in fv|>typeof|>fieldnames)...)
Base.reverse(fv::CompositeVector)=typeof(fv)((name==:contents ? reverse(fv.contents) : getfield(fv,name) for name in fv|>typeof|>fieldnames)...)
function Base.similar(fv::CompositeVector,dtype::Type{T}=eltype(fv),dims::NTuple{N,Int}=size(fv)) where {T,N}
    typeof(fv).name.wrapper((name==:contents ? similar(fv.contents,dtype,dims) : getfield(fv,name) for name in fv|>typeof|>fieldnames)...)
end
Base.iterate(fv::CompositeVector,state=1)=iterate(fv.contents,state)
Base.keys(fv::CompositeVector)=keys(fv.contents)
Base.values(fv::CompositeVector)=values(fv.contents)
Base.pairs(fv::CompositeVector)=pairs(fv.contents)

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
Base.:(==)(cd1::CD,cd2::CD) where CD<:CompositeDict=all(getfield(cd1,name)==getfield(cd2,name) for name in CD|>fieldnames)
Base.isequal(cd1::CD,cd2::CD) where CD<:CompositeDict=all(isequal(getfield(cd1,name),getfield(cd2,name)) for name in CD|>fieldnames)
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

end #module
