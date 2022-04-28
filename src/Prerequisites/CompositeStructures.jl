module CompositeStructures

using ..Traits: getcontent, dissolve, rawtype, efficientoperations

import ..Traits: contentnames, dissolve

export CompositeTuple, CompositeNTuple, CompositeVector, CompositeDict
export NamedContainer

"""
    CompositeTuple{T<:Tuple}

A composite tuple can be considered as a tuple that is implemented by including an ordinary `Tuple` as its data attribute.
"""
abstract type CompositeTuple{T<:Tuple} end
@inline contentnames(::Type{<:CompositeTuple}) = (:contents,)
@inline dissolve(ct::CompositeTuple, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(ct, :contents), args...; kwargs...)

@inline Base.length(::CompositeTuple{T}) where T = fieldcount(T)
@inline Base.length(::Type{<:CompositeTuple{T}}) where T = fieldcount(T)
@inline Base.eltype(::CompositeTuple{T}) where T = eltype(T)
@inline Base.eltype(::Type{<:CompositeTuple{T}}) where T = eltype(T)
@inline Base.hash(ct::CompositeTuple, h::UInt) = hash(dissolve(ct), h)
@inline Base.:(==)(ct1::CompositeTuple, ct2::CompositeTuple) = ==(efficientoperations, ct1, ct2)
@inline Base.isequal(ct1::CompositeTuple, ct2::CompositeTuple) = isequal(efficientoperations, ct1, ct2)
@inline Base.getindex(ct::CompositeTuple, i::Union{<:Integer, CartesianIndex}) = getcontent(ct, :contents)[i]
@inline Base.getindex(ct::CompositeTuple, inds) = rawtype(typeof(ct))(dissolve(ct, getindex, (inds,))...)
@inline Base.lastindex(ct::CompositeTuple) = lastindex(getcontent(ct, :contents))
@inline Base.iterate(ct::CompositeTuple) = iterate(getcontent(ct, :contents))
@inline Base.iterate(ct::CompositeTuple, state) = iterate(getcontent(ct, :contents), state)
@inline Base.iterate(rv::Iterators.Reverse{<:CompositeTuple}, state=length(rv.itr)) = (state < 1) ? nothing : (rv.itr[state], state-1)
@inline Base.keys(ct::CompositeTuple) = keys(getcontent(ct, :contents))
@inline Base.values(ct::CompositeTuple) = values(getcontent(ct, :contents))
@inline Base.pairs(ct::CompositeTuple) = pairs(getcontent(ct, :contents))
@inline Base.reverse(ct::CompositeTuple) = rawtype(typeof(ct))(dissolve(ct, reverse)...)
@inline Base.Tuple(ct::CompositeTuple) = getcontent(ct, :contents)

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
@inline contentnames(::Type{<:CompositeVector}) = (:contents,)
@inline dissolve(cv::CompositeVector, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(cv, :contents), args...; kwargs...)

@inline Base.size(cv::CompositeVector) = size(getcontent(cv, :contents))
@inline Base.size(cv::CompositeVector, i) = size(getcontent(cv, :contents), i)
@inline Base.length(cv::CompositeVector) = length(getcontent(cv, :contents))
@inline Base.:(==)(cv1::CompositeVector, cv2::CompositeVector) = ==(efficientoperations, cv1, cv2)
@inline Base.isequal(cv1::CompositeVector, cv2::CompositeVector) = isequal(efficientoperations, cv1, cv2)
@inline Base.getindex(cv::CompositeVector, i::Union{<:Integer, CartesianIndex}) = getcontent(cv, :contents)[i]
@inline Base.getindex(cv::CompositeVector, inds) = rawtype(typeof(cv))(dissolve(cv, getindex, (inds,))...)
@inline Base.lastindex(cv::CompositeVector) = lastindex(getcontent(cv, :contents))
@inline Base.setindex!(cv::CompositeVector, value, inds) = (getcontent(cv, :contents)[inds] = value)
@inline Base.push!(cv::CompositeVector, values...) = (push!(getcontent(cv, :contents), values...); cv)
@inline Base.pushfirst!(cv::CompositeVector, values...) = (pushfirst!(getcontent(cv, :contents), values...); cv)
@inline Base.insert!(cv::CompositeVector, index::Integer, value) = (insert!(getcontent(cv, :contents), index, value); cv)
@inline Base.append!(cv::CompositeVector, values) = (append!(getcontent(cv, :contents), values); cv)
@inline Base.prepend!(cv::CompositeVector, values) = (prepend!(getcontent(cv, :contents), values); cv)
@inline Base.splice!(cv::CompositeVector, index::Integer, replacement=Base._default_splice) = splice!(getcontent(cv, :contents), index, replacement)
@inline Base.splice!(cv::CompositeVector, range::UnitRange{<:Integer}, replacement=Base._default_splice) = rawtype(typeof(cv))(dissolve(cv, splice!, (range, replacement))...)
@inline Base.deleteat!(cv::CompositeVector, indices) = (deleteat!(getcontent(cv, :contents), indices); cv)
@inline Base.pop!(cv::CompositeVector) = pop!(getcontent(cv, :contents))
@inline Base.popfirst!(cv::CompositeVector) = popfirst!(getcontent(cv, :contents))
@inline Base.empty!(cv::CompositeVector) = (empty!(getcontent(cv, :contents)); cv)
@inline Base.empty(cv::CompositeVector) = rawtype(typeof(cv))(dissolve(cv, empty)...)
@inline Base.reverse(cv::CompositeVector) =  rawtype(typeof(cv))(dissolve(cv, reverse)...)
@inline Base.similar(cv::CompositeVector, dtype::Type=eltype(cv), dims::Tuple{Vararg{Int}}=size(cv)) = rawtype(typeof(cv))(dissolve(cv, similar, (dtype, dims))...)
@inline Base.iterate(cv::CompositeVector, state=1) = iterate(getcontent(cv, :contents), state)
@inline Base.keys(cv::CompositeVector) = keys(getcontent(cv, :contents))
@inline Base.values(cv::CompositeVector) = values(getcontent(cv, :contents))
@inline Base.pairs(cv::CompositeVector) = pairs(getcontent(cv, :contents))

"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
@inline contentnames(::Type{<:CompositeDict}) = (:contents,)
@inline dissolve(cd::CompositeDict, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(cd, :contents), args...; kwargs...)

@inline Base.isempty(cd::CompositeDict) = isempty(getcontent(cd, :contents))
@inline Base.length(cd::CompositeDict) = length(getcontent(cd, :contents))
@inline Base.haskey(cd::CompositeDict, key) = haskey(getcontent(cd, :contents), key)
@inline Base.in(p::Pair, cd::CompositeDict, valcmp=(==)) = in(p, getcontent(cd, :contents), valcmp)
@inline Base.:(==)(cd1::CompositeDict, cd2::CompositeDict) = ==(efficientoperations, cd1, cd2)
@inline Base.isequal(cd1::CompositeDict, cd2::CompositeDict) = isequal(efficientoperations, cd1, cd2)
@inline Base.get(cd::CompositeDict, key, default) = get(getcontent(cd, :contents), key, default)
@inline Base.get(f::Base.Callable, cd::CompositeDict, key) = get(f, getcontent(cd, :contents), key)
@inline Base.getkey(cd::CompositeDict, key, default) = getkey(getcontent(cd, :contents), key, default)
@inline Base.getindex(cd::CompositeDict{K, V}, index::K) where {K, V} = getcontent(cd, :contents)[index]
@inline Base.push!(cd::CompositeDict, ps::Pair...) = (push!(getcontent(cd, :contents), ps...); cd)
@inline Base.get!(cd::CompositeDict, key, default) = get!(getcontent(cd, :contents), key, default)
@inline Base.get!(default::Union{Function, Type}, cd::CompositeDict, key) = get!(default, getcontent(cd, :contents), key)
@inline Base.setindex!(cd::CompositeDict{K, V}, value::V, index::K) where {K, V} = (getcontent(cd, :contents)[index] = value)
@inline Base.pop!(cd::CompositeDict) = pop!(getcontent(cd, :contents))
@inline Base.pop!(cd::CompositeDict, key) = pop!(getcontent(cd, :contents), key)
@inline Base.pop!(cd::CompositeDict, key, default) = pop!(getcontent(cd, :contents), key, default)
@inline Base.delete!(cd::CompositeDict, key) = (delete!(getcontent(cd, :contents), key); cd)
@inline Base.empty!(cd::CompositeDict) = (empty!(getcontent(cd, :contents)); cd)
@inline Base.merge(cd::CD, others::CD...) where CD <: CompositeDict = merge!(empty(cd), cd, others...)
@inline Base.merge(combine::Function, cd::CD, others::CD...) where {CD<:CompositeDict} = merge!(combine, empty(cd), cd, others...)
@inline Base.empty(cd::CompositeDict) = rawtype(typeof(cd))(dissolve(cd, empty)...)
@inline Base.iterate(cd::CompositeDict) = iterate(getcontent(cd, :contents))
@inline Base.iterate(cd::CompositeDict, state) = iterate(getcontent(cd, :contents), state)
@inline Base.keys(cd::CompositeDict) = keys(getcontent(cd, :contents))
@inline Base.values(cd::CompositeDict) = values(getcontent(cd, :contents))
@inline Base.pairs(cd::CompositeDict) = pairs(getcontent(cd, :contents))

"""
    NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}

NamedContainer is just a wrapper of Julia NamedTuple, but not a composite type.
"""
const NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}
@inline NamedContainer{Names}(contents::Tuple) where Names = NamedTuple{Names}(contents)

end #module
