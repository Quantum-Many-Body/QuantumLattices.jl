module CompositeStructures

using ..Traits: getcontent, disolve, rawtype, efficientoperations

import ..Traits: contentnames 

export CompositeTuple, CompositeNTuple, CompositeVector, CompositeDict
export NamedContainer

"""
    CompositeTuple{T<:Tuple}

A composite tuple can be considered as a tuple that is implemented by including an ordinary `Tuple` as its data attribute.
"""
abstract type CompositeTuple{T<:Tuple} end
contentnames(::Type{<:CompositeTuple}) = (:contents,)

Base.length(::CompositeTuple{T}) where T = fieldcount(T)
Base.length(::Type{<:CompositeTuple{T}}) where T = fieldcount(T)
Base.eltype(::CompositeTuple{T}) where T = eltype(T)
Base.eltype(::Type{<:CompositeTuple{T}}) where T = eltype(T)
Base.hash(ct::CompositeTuple, h::UInt) = hash(disolve(ct), h)
Base.:(==)(ct1::CompositeTuple, ct2::CompositeTuple) = ==(efficientoperations, ct1, ct2)
Base.isequal(ct1::CompositeTuple, ct2::CompositeTuple) = isequal(efficientoperations, ct1, ct2)
Base.getindex(ct::CompositeTuple, i::Union{<:Integer, CartesianIndex}) = getcontent(ct, :contents)[i]
Base.getindex(ct::CompositeTuple, inds) = rawtype(typeof(ct))(disolve(ct, getindex, (inds,))...)
Base.lastindex(ct::CompositeTuple) = lastindex(getcontent(ct, :contents))
Base.iterate(ct::CompositeTuple) = iterate(getcontent(ct, :contents))
Base.iterate(ct::CompositeTuple, state) = iterate(getcontent(ct, :contents), state)
Base.iterate(rv::Iterators.Reverse{<:CompositeTuple}, state=length(rv.itr)) = (state < 1) ? nothing : (rv.itr[state], state-1)
Base.keys(ct::CompositeTuple) = keys(getcontent(ct, :contents))
Base.values(ct::CompositeTuple) = values(getcontent(ct, :contents))
Base.pairs(ct::CompositeTuple) = pairs(getcontent(ct, :contents))
Base.reverse(ct::CompositeTuple) = rawtype(typeof(ct))(disolve(ct, reverse)...)
Base.convert(::Type{Tuple}, ct::CompositeTuple) = getcontent(ct, :contents)

"""
    CompositeNTuple{N, T}

A composite ntuple can be considered as a ntuple that is implemented by including an ordinary `NTuple` as its data attribute.

Alias for `CompositeTuple{NTuple{N, T}}`.
"""
const CompositeNTuple{N, T} = CompositeTuple{NTuple{N, T}}

"""
    CompositeVector{T}

A composite vector can be considered as a vector that is implemented by including a concrete subtype of `AbstractVector` as its data attribute.
"""
abstract type CompositeVector{T} <: AbstractVector{T} end
contentnames(::Type{<:CompositeVector}) = (:contents,)

Base.size(cv::CompositeVector) = size(getcontent(cv, :contents))
Base.size(cv::CompositeVector, i) = size(getcontent(cv, :contents), i)
Base.length(cv::CompositeVector) = length(getcontent(cv, :contents))
Base.:(==)(cv1::CompositeVector, cv2::CompositeVector) = ==(efficientoperations, cv1, cv2)
Base.isequal(cv1::CompositeVector, cv2::CompositeVector) = isequal(efficientoperations, cv1, cv2)
Base.getindex(cv::CompositeVector, i::Union{<:Integer, CartesianIndex}) = getcontent(cv, :contents)[i]
Base.getindex(cv::CompositeVector, inds) = rawtype(typeof(cv))(disolve(cv, getindex, (inds,))...)
Base.lastindex(cv::CompositeVector) = lastindex(getcontent(cv, :contents))
Base.setindex!(cv::CompositeVector, value, inds) = (getcontent(cv, :contents)[inds] = value)
Base.push!(cv::CompositeVector, values...) = (push!(getcontent(cv, :contents), values...); cv)
Base.pushfirst!(cv::CompositeVector, values...) = (pushfirst!(getcontent(cv, :contents), values...); cv)
Base.insert!(cv::CompositeVector, index::Integer, value) = (insert!(getcontent(cv, :contents), index, value); cv)
Base.append!(cv::CompositeVector, values) = (append!(getcontent(cv, :contents), values); cv)
Base.prepend!(cv::CompositeVector, values) = (prepend!(getcontent(cv, :contents), values); cv)
Base.splice!(cv::CompositeVector, index::Integer, replacement=Base._default_splice) = splice!(getcontent(cv, :contents), index, replacement)
Base.splice!(cv::CompositeVector, range::UnitRange{<:Integer}, replacement=Base._default_splice) = rawtype(typeof(cv))(disolve(cv, splice!, (range, replacement))...)
Base.deleteat!(cv::CompositeVector, indices) = (deleteat!(getcontent(cv, :contents), indices); cv)
Base.pop!(cv::CompositeVector) = pop!(getcontent(cv, :contents))
Base.popfirst!(cv::CompositeVector) = popfirst!(getcontent(cv, :contents))
Base.empty!(cv::CompositeVector) = (empty!(getcontent(cv, :contents)); cv)
Base.empty(cv::CompositeVector) = rawtype(typeof(cv))(disolve(cv, empty)...)
Base.reverse(cv::CompositeVector) =  rawtype(typeof(cv))(disolve(cv, reverse)...)
Base.similar(cv::CompositeVector, dtype::Type=eltype(cv), dims::Tuple{Vararg{Int}}=size(cv)) = rawtype(typeof(cv))(disolve(cv, similar, (dtype, dims))...)
Base.iterate(cv::CompositeVector, state=1) = iterate(getcontent(cv, :contents), state)
Base.keys(cv::CompositeVector) = keys(getcontent(cv, :contents))
Base.values(cv::CompositeVector) = values(getcontent(cv, :contents))
Base.pairs(cv::CompositeVector) = pairs(getcontent(cv, :contents))
Base.convert(::Type{Vector}, cv::CompositeVector) = getcontent(cv, :contents)

"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
contentnames(::Type{<:CompositeDict}) = (:contents,)

Base.isempty(cd::CompositeDict) = isempty(getcontent(cd, :contents))
Base.length(cd::CompositeDict) = length(getcontent(cd, :contents))
Base.haskey(cd::CompositeDict, key) = haskey(getcontent(cd, :contents), key)
Base.in(p::Pair, cd::CompositeDict, valcmp=(==)) = in(p, getcontent(cd, :contents), valcmp)
Base.:(==)(cd1::CompositeDict, cd2::CompositeDict) = ==(efficientoperations, cd1, cd2)
Base.isequal(cd1::CompositeDict, cd2::CompositeDict) = isequal(efficientoperations, cd1, cd2)
Base.get(cd::CompositeDict, key, default) = get(getcontent(cd, :contents), key, default)
Base.get(f::Base.Callable, cd::CompositeDict, key) = get(f, getcontent(cd, :contents), key)
Base.getkey(cd::CompositeDict, key, default) = getkey(getcontent(cd, :contents), key, default)
Base.getindex(cd::CompositeDict{K, V}, index::K) where {K, V} = getcontent(cd, :contents)[index]
Base.push!(cd::CompositeDict, ps::Pair...) = (push!(getcontent(cd, :contents), ps...); cd)
Base.get!(cd::CompositeDict, key, default) = get!(getcontent(cd, :contents), key, default)
Base.get!(default::Union{Function, Type}, cd::CompositeDict, key) = get!(default, getcontent(cd, :contents), key)
Base.setindex!(cd::CompositeDict{K, V}, value::V, index::K) where {K, V} = (getcontent(cd, :contents)[index] = value)
Base.pop!(cd::CompositeDict) = pop!(getcontent(cd, :contents))
Base.pop!(cd::CompositeDict, key) = pop!(getcontent(cd, :contents), key)
Base.pop!(cd::CompositeDict, key, default) = pop!(getcontent(cd, :contents), key, default)
Base.delete!(cd::CompositeDict, key) = (delete!(getcontent(cd, :contents), key); cd)
Base.empty!(cd::CompositeDict) = (empty!(getcontent(cd, :contents)); cd)
Base.merge(cd::CD, others::CD...) where CD <: CompositeDict = merge!(empty(cd), cd, others...)
Base.merge(combine::Function, cd::CD, others::CD...) where {CD<:CompositeDict} = merge!(combine, empty(cd), cd, others...)
Base.empty(cd::CompositeDict) = rawtype(typeof(cd))(disolve(cd, empty)...)
Base.iterate(cd::CompositeDict) = iterate(getcontent(cd, :contents))
Base.iterate(cd::CompositeDict, state) = iterate(getcontent(cd, :contents), state)
Base.keys(cd::CompositeDict) = keys(getcontent(cd, :contents))
Base.values(cd::CompositeDict) = values(getcontent(cd, :contents))
Base.pairs(cd::CompositeDict) = pairs(getcontent(cd, :contents))
Base.convert(::Type{Dict}, cd::CompositeDict) = getcontent(cd, :contents)

"""
    NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}

NamedContainer is just a wrapper of Julia NamedTuple, but not a composite type.
"""
const NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}

"""
    NamedContainer{Names}(contents) where Names -> NamedTuple{Names, typeof(contents)}

Construct a named container.
"""
@generated function NamedContainer{Names}(contents::Tuple) where Names
    @assert length(Names) == fieldcount(contents) "NamedContainer error: dismatched length between names and contents."
    fieldcount(contents) == 0 && return NamedTuple()
    return Expr(:tuple, [:($name = contents[$i]) for (i, name) in enumerate(Names)]...)
end

end #module
