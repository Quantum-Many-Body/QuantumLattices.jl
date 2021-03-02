module CompositeStructures

using ..TypeTraits: efficientoperations, rawtype, Field

export CompositeTuple, CompositeNTuple, CompositeVector, CompositeDict
export NamedContainer

"""
    CompositeTuple{T<:Tuple}

A composite tuple can be considered as a tuple that is implemented by including an ordinary `Tuple` as its data attribute.
"""
abstract type CompositeTuple{T<:Tuple} end
Base.fieldnames(::Type{Field}, ::Type{<:CompositeTuple}) = (:contents,)
Base.map(f::Function, ct::CompositeTuple, a::Field{:contents}, args::Tuple, kwargs::NamedTuple) = f(getproperty(ct, a), args...; kwargs...)

Base.length(::CompositeTuple{T}) where T = fieldcount(T)
Base.length(::Type{<:CompositeTuple{T}}) where T = fieldcount(T)
Base.eltype(::CompositeTuple{T}) where T = eltype(T)
Base.eltype(::Type{<:CompositeTuple{T}}) where T = eltype(T)
Base.hash(ct::CompositeTuple, h::UInt) = hash(Tuple(Field, ct), h)
Base.:(==)(ct1::CompositeTuple, ct2::CompositeTuple) = ==(efficientoperations, ct1, ct2)
Base.isequal(ct1::CompositeTuple, ct2::CompositeTuple) = isequal(efficientoperations, ct1, ct2)
Base.getindex(ct::CompositeTuple, i::Union{<:Integer, CartesianIndex}) = getproperty(ct, Field(:contents))[i]
Base.getindex(ct::CompositeTuple, inds) = rawtype(typeof(ct))(Tuple(Field, ct, getindex, (inds,))...)
Base.lastindex(ct::CompositeTuple) = lastindex(getproperty(ct, Field(:contents)))
Base.iterate(ct::CompositeTuple) = iterate(getproperty(ct, Field(:contents)))
Base.iterate(ct::CompositeTuple, state) = iterate(getproperty(ct, Field(:contents)), state)
Base.iterate(rv::Iterators.Reverse{<:CompositeTuple}, state=length(rv.itr)) = (state < 1) ? nothing : (rv.itr[state], state-1)
Base.keys(ct::CompositeTuple) = keys(getproperty(ct, Field(:contents)))
Base.values(ct::CompositeTuple) = values(getproperty(ct, Field(:contents)))
Base.pairs(ct::CompositeTuple) = pairs(getproperty(ct, Field(:contents)))
Base.reverse(ct::CompositeTuple) = rawtype(typeof(ct))(Tuple(Field, ct, reverse)...)
Base.convert(::Type{Tuple}, ct::CompositeTuple) = getproperty(ct, Field(:contents))

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
Base.fieldnames(::Type{Field}, ::Type{<:CompositeVector}) = (:contents,)
Base.map(f::Function, cv::CompositeVector, a::Field{:contents}, args::Tuple, kwargs::NamedTuple) = f(getproperty(cv, a), args...; kwargs...)

Base.size(cv::CompositeVector) = size(getproperty(cv, Field(:contents)))
Base.size(cv::CompositeVector, i) = size(getproperty(cv, Field(:contents)), i)
Base.length(cv::CompositeVector) = length(getproperty(cv, Field(:contents)))
Base.:(==)(cv1::CompositeVector, cv2::CompositeVector) = ==(efficientoperations, cv1, cv2)
Base.isequal(cv1::CompositeVector, cv2::CompositeVector) = isequal(efficientoperations, cv1, cv2)
Base.getindex(cv::CompositeVector, i::Union{<:Integer, CartesianIndex}) = getproperty(cv, Field(:contents))[i]
Base.getindex(cv::CompositeVector, inds) = rawtype(typeof(cv))(Tuple(Field, cv, getindex, (inds,))...)
Base.lastindex(cv::CompositeVector) = lastindex(getproperty(cv, Field(:contents)))
Base.setindex!(cv::CompositeVector, value, inds) = (getproperty(cv, Field(:contents))[inds] = value)
Base.push!(cv::CompositeVector, values...) = (push!(getproperty(cv, Field(:contents)), values...); cv)
Base.pushfirst!(cv::CompositeVector, values...) = (pushfirst!(getproperty(cv, Field(:contents)), values...); cv)
Base.insert!(cv::CompositeVector, index::Integer, value) = (insert!(getproperty(cv, Field(:contents)), index, value); cv)
Base.append!(cv::CompositeVector, values) = (append!(getproperty(cv, Field(:contents)), values); cv)
Base.prepend!(cv::CompositeVector, values) = (prepend!(getproperty(cv, Field(:contents)), values); cv)
Base.splice!(cv::CompositeVector, index::Integer, replacement=Base._default_splice) = splice!(getproperty(cv, Field(:contents)), index, replacement)
Base.splice!(cv::CompositeVector, range::UnitRange{<:Integer}, replacement=Base._default_splice) = rawtype(typeof(cv))(Tuple(Field, cv, splice!, (range, replacement))...)
Base.deleteat!(cv::CompositeVector, indices) = (deleteat!(getproperty(cv, Field(:contents)), indices); cv)
Base.pop!(cv::CompositeVector) = pop!(getproperty(cv, Field(:contents)))
Base.popfirst!(cv::CompositeVector) = popfirst!(getproperty(cv, Field(:contents)))
Base.empty!(cv::CompositeVector) = (empty!(getproperty(cv, Field(:contents))); cv)
Base.empty(cv::CompositeVector) = rawtype(typeof(cv))(Tuple(Field, cv, empty)...)
Base.reverse(cv::CompositeVector) =  rawtype(typeof(cv))(Tuple(Field, cv, reverse)...)
Base.similar(cv::CompositeVector, dtype::Type=eltype(cv), dims::Tuple{Vararg{Int}}=size(cv)) = rawtype(typeof(cv))(Tuple(Field, cv, similar, (dtype, dims))...)
Base.iterate(cv::CompositeVector, state=1) = iterate(getproperty(cv, Field(:contents)), state)
Base.keys(cv::CompositeVector) = keys(getproperty(cv, Field(:contents)))
Base.values(cv::CompositeVector) = values(getproperty(cv, Field(:contents)))
Base.pairs(cv::CompositeVector) = pairs(getproperty(cv, Field(:contents)))
Base.convert(::Type{Vector}, cv::CompositeVector) = getproperty(cv, Field(:contents))

"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
Base.fieldnames(::Type{Field}, ::Type{<:CompositeDict}) = (:contents,)
Base.map(f::Function, cv::CompositeDict, a::Field{:contents}, args::Tuple, kwargs::NamedTuple) = f(getproperty(cv, a), args...; kwargs...)

Base.isempty(cd::CompositeDict) = isempty(getproperty(cd, Field(:contents)))
Base.length(cd::CompositeDict) = length(getproperty(cd, Field(:contents)))
Base.haskey(cd::CompositeDict, key) = haskey(getproperty(cd, Field(:contents)), key)
Base.in(p::Pair, cd::CompositeDict, valcmp=(==)) = in(p, getproperty(cd, Field(:contents)), valcmp)
Base.:(==)(cd1::CompositeDict, cd2::CompositeDict) = ==(efficientoperations, cd1, cd2)
Base.isequal(cd1::CompositeDict, cd2::CompositeDict) = isequal(efficientoperations, cd1, cd2)
Base.get(cd::CompositeDict, key, default) = get(getproperty(cd, Field(:contents)), key, default)
Base.get(f::Base.Callable, cd::CompositeDict, key) = get(f, getproperty(cd, Field(:contents)), key)
Base.getkey(cd::CompositeDict, key, default) = getkey(getproperty(cd, Field(:contents)), key, default)
Base.getindex(cd::CompositeDict{K, V}, index::K) where {K, V} = getproperty(cd, Field(:contents))[index]
Base.push!(cd::CompositeDict, ps::Pair...) = (push!(getproperty(cd, Field(:contents)), ps...); cd)
Base.get!(cd::CompositeDict, key, default) = get!(getproperty(cd, Field(:contents)), key, default)
Base.get!(default::Union{Function, Type}, cd::CompositeDict, key) = get!(default, getproperty(cd, Field(:contents)), key)
Base.setindex!(cd::CompositeDict{K, V}, value::V, index::K) where {K, V} = (getproperty(cd, Field(:contents))[index] = value)
Base.pop!(cd::CompositeDict) = pop!(getproperty(cd, Field(:contents)))
Base.pop!(cd::CompositeDict, key) = pop!(getproperty(cd, Field(:contents)), key)
Base.pop!(cd::CompositeDict, key, default) = pop!(getproperty(cd, Field(:contents)), key, default)
Base.delete!(cd::CompositeDict, key) = (delete!(getproperty(cd, Field(:contents)), key); cd)
Base.empty!(cd::CompositeDict) = (empty!(getproperty(cd, Field(:contents))); cd)
Base.merge(cd::CD, others::CD...) where CD <: CompositeDict = merge!(empty(cd), cd, others...)
Base.merge(combine::Function, cd::CD, others::CD...) where {CD<:CompositeDict} = merge!(combine, empty(cd), cd, others...)
Base.empty(cd::CompositeDict) = rawtype(typeof(cd))(Tuple(Field, cd, empty)...)
Base.iterate(cd::CompositeDict) = iterate(getproperty(cd, Field(:contents)))
Base.iterate(cd::CompositeDict, state) = iterate(getproperty(cd, Field(:contents)), state)
Base.keys(cd::CompositeDict) = keys(getproperty(cd, Field(:contents)))
Base.values(cd::CompositeDict) = values(getproperty(cd, Field(:contents)))
Base.pairs(cd::CompositeDict) = pairs(getproperty(cd, Field(:contents)))
Base.convert(::Type{Dict}, cd::CompositeDict) = getproperty(cd, Field(:contents))

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
