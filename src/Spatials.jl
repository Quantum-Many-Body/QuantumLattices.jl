module Spatials

using Base.Iterators: flatten, product
using DelimitedFiles: writedlm
using LinearAlgebra: cross, det, dot, norm
using NearestNeighbors: KDTree, inrange, knn
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series, @layout
using StaticArrays: MVector, SVector
using ..QuantumLattices: OneAtLeast, OneOrMore, rank
using ..QuantumNumbers: Momenta, 𝕂, period, periods
using ..Toolkit: atol, rtol, efficientoperations, CompositeDict, DirectProductedVectorSpace, Float, Segment, VectorSpace, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceEnumerative, VectorSpaceStyle, concatenate, getcontent, subscript

import ..QuantumLattices: decompose, dimension, expand, kind, matrix, shape
import ..QuantumNumbers: 𝕂¹, 𝕂², 𝕂³
import ..QuantumOperators: scalartype
import ..Toolkit: contentnames

export azimuth, azimuthd, direction, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, FractionalReciprocalSpace, Lattice, Neighbors, Point, ProductedReciprocalSpace, ReciprocalCurve, ReciprocalScatter, ReciprocalSpace, ReciprocalZone, ReciprocalPath
export bonds, bonds!, fractionals, icoordinate, isintracell, label, nneighbor, rcoordinate, savedlm, selectpath, shrink, ticks
export hexagon120°map, hexagon60°map, linemap, rectanglemap, @hexagon_str, @line_str, @rectangle_str

"""
    distance(p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}) -> Number

Get the distance between two points.

!!! note
    Compared to `norm(p₁-p₂)`, this function avoids the memory allocation for `p₁-p₂`, thus is more efficient.
"""
function distance(p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number})
    result = zero(promote_type(eltype(p₁), eltype(p₂)))
    for i in eachindex(p₁, p₂)
        result = result + (p₁[i]-p₂[i])^2
    end
    return sqrt(result)
end

"""
    azimuth(v::AbstractVector{<:Number}) -> Number

Get the azimuth angle in radians of a vector.
"""
function azimuth(v::AbstractVector{<:Number})
    @assert length(v)∈(1, 2, 3) "azimuth error: wrong dimensioned input vector."
    if length(v)==1
        result = acos(first(v)/abs(first(v)))
    else
        e₁, e₂ = v
        result = acos(e₁/(length(v)==3 ? sqrt(e₁^2+e₂^2) : norm(v)))
        e₂<0 && (result = 2*convert(typeof(result), pi) - result)
    end
    isnan(result) && (result = zero(result))
    return result
end

"""
    azimuthd(v::AbstractVector{<:Number}) -> Number

Get the azimuth angle in degrees of a vector.
"""
function azimuthd(v::AbstractVector{<:Number})
    @assert length(v)∈(1, 2, 3) "azimuthd error: wrong dimensioned input vector."
    if length(v)==1
        result = acosd(first(v)/abs(first(v)))
    else
        e₁, e₂ = v
        result = acosd(e₁/(length(v)==3 ? sqrt(e₁^2+e₂^2) : norm(v)))
        e₂<0 && (result = 360 - result)
    end
    isnan(result) && (result = zero(result))
    return result
end

"""
    polar(v::AbstractVector{<:Number}) -> Number

Get the polar angle in radians of a vector.
"""
function polar(v::AbstractVector{<:Number})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acos(v[end]/norm(v))
end

"""
    polard(v::AbstractVector{<:Number}) -> Number

Get the polar angle in degrees of a vector.
"""
function polard(v::AbstractVector{<:Number})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acosd(v[end]/norm(v))
end

"""
    direction(v::Char, args...) -> SVector{3, <:Number}
    direction(v::Number, unit::Symbol) -> SVector{2, <:Number}
    direction(v::Tuple{Number, Number}, unit::Symbol) -> SVector{3, <:Number}
    direction(v::AbstractVector{<:Number}, args...) -> AbstractVector{<:Number}

Get the unit vector that specifies the direction of a vector.
"""
function direction(v::Char, args...)
    @assert lowercase(v)∈('x', 'y', 'z') "direction error: not one of 'x', 'y' and 'z'."
    return lowercase(v)=='x' ? SVector(1, 0, 0) : lowercase(v)=='y' ? SVector(0, 1, 0) : SVector(0, 0, 1)
end
function direction(v::Number, unit::Symbol)
    @assert unit∈(:degree, :radian) "direction error: not supported unit of angles."
    if unit==:degree
        return SVector(cosd(v), sind(v))
    else
        return SVector(cos(v), sin(v))
    end
end
function direction(v::Tuple{Number, Number}, unit::Symbol)
    @assert unit∈(:degree, :radian) "direction error: not supported unit of angles."
    θ, ϕ = v[1], v[2]
    if unit==:degree
        sinθ = sind(θ)
        return SVector(sinθ*cosd(ϕ), sinθ*sind(ϕ), cosd(θ))
    else
        sinθ = sin(θ)
        return SVector(sinθ*cos(ϕ), sinθ*sin(ϕ), cos(θ))
    end
end
function direction(v::AbstractVector{<:Number}, args...)
    n = norm(v)
    return map(i->i/n, v)
end

"""
    volume(vectors::AbstractVector{<:AbstractVector{<:Number}}) -> Number
    volume(v::AbstractVector{<:Number}) -> Number
    volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}) -> Number
    volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number}) -> Number

Get the volume spanned by the input vectors.
"""
function volume(vectors::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(vectors)∈(1, 2, 3) "volume error: unsupported number of vectors."
    length(vectors)==1 && return volume(first(vectors))
    length(vectors)==2 && return volume(first(vectors), last(vectors))
    v₁, v₂, v₃ = vectors
    return volume(v₁, v₂, v₃)
end
@inline volume(v::AbstractVector{<:Number}) = norm(v)
function volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number})
    @assert length(v₁)==length(v₂) "volume error: mismatched dimension of vectors."
    @assert length(v₁)∈(1, 2, 3) "volume error: unsupported dimension of vectors."
    length(v₁)==1 && return zero(eltype(v₁))
    length(v₁)==2 && return abs(det(hcat(SVector{2}(v₁), SVector{2}(v₂))))
    return norm(cross(SVector{3}(v₁), SVector{3}(v₂)))
end
function volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number})
    @assert length(v₁)==length(v₂)==length(v₃) "volume error: mismatched dimension of vectors."
    @assert length(v₁)∈(1, 2, 3) "volume error: unsupported dimension of vectors."
    length(v₁)∈(1, 2) && return zero(eltype(v₁))
    return abs(det(hcat(SVector{3}(v₁), SVector{3}(v₂), SVector{3}(v₃))))
end

"""
    decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}) -> Number
    decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}) -> Tuple{Number, Number}
    decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number}) -> Tuple{Number, Number, Number}
    decompose(v₀::AbstractVector{<:Number}, vs::AbstractVector{<:AbstractVector{<:Number}}) -> Vector{<:Number}

Decompose a vector with respect to input basis vectors.
"""
function decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number})
    @assert length(v₀)==length(v₁) "decompose error: mismatched length of input vectors."
    n₀, n₁ = norm(v₀), norm(v₁)
    abs(n₀)≈0 && return zero(n₀)
    sign = dot(v₀, v₁) / n₀ / n₁
    @assert isapprox(abs(sign), 1.0, atol=atol, rtol=rtol) "decompose error: insufficient basis vectors."
    return sign*n₀/n₁
end
function decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number})
    @assert length(v₀)==length(v₁)==length(v₂) "decompose error: mismatched length of input vectors."
    @assert 1<length(v₀)<4 "decompose error: unsupported dimension($(length(v₀))) of input vectors."
    if length(v₀) == 2
        result = inv(hcat(SVector{2}(v₁), SVector{2}(v₂))) * SVector{2}(v₀)
        return result[1], result[2]
    else
        v₁, v₂ = SVector{3}(v₁), SVector{3}(v₂)
        x₁, x₂, x₃ = decompose(v₀, v₁, v₂, cross(v₁, v₂))
        @assert isapprox(x₃, 0.0, atol=atol, rtol=rtol) "decompose error: insufficient basis vectors."
        return x₁, x₂
    end
end
function decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number})
    @assert length(v₀)==length(v₁)==length(v₂)==length(v₃) "decompose error: mismatched length of input vectors."
    @assert length(v₀)==3 "decompose error: unsupported dimension($(length(v₀))) of input vectors."
    result = inv(hcat(SVector{3}(v₁), SVector{3}(v₂), SVector{3}(v₃))) * SVector{3}(v₀)
    return result[1], result[2], result[3]
end
function decompose(v₀::AbstractVector{<:Number}, vs::AbstractVector{<:AbstractVector{<:Number}})
    @assert 0<length(vs)<4 "decompose error: not supported number ($(length(vs))) of vectors."
    if length(vs)==1
        result = [decompose(v₀, first(vs))]
    elseif length(vs)==2
        result = collect(decompose(v₀, first(vs), last(vs)))
    else
        v₁, v₂, v₃ = vs
        result = collect(decompose(v₀, v₁, v₂, v₃))
    end
    return result
end

"""
    decompose(m::AbstractMatrix{<:Number}, m₀::AbstractMatrix{<:Number}) -> Number
    decompose(m::AbstractMatrix{<:Number}, ms::Tuple{Vararg{AbstractMatrix{<:Number}}}) -> Tuple{Vararg{Number}}
    decompose(m::AbstractMatrix{<:Number}, ms::AbstractVector{<:AbstractMatrix{<:Number}}) -> Vector{<:Number}

Decompose a matrix.
"""
function decompose(m::AbstractMatrix{<:Number}, m₀::AbstractMatrix{<:Number})
    @assert size(m)==size(m₀) "decompose error: mismatched size."
    result = zero(promote_type(eltype(m), eltype(m₀)))
    n = zero(result)
    for (i, j) in zip(m, m₀)
        result += conj(j) * i
        n += conj(j) * j
    end
    return result/n
end
@inline decompose(m::AbstractMatrix{<:Number}, ms::Tuple{Vararg{AbstractMatrix{<:Number}}}) = map(m₀->decompose(m, m₀), ms)
@inline decompose(m::AbstractMatrix{<:Number}, ms::AbstractVector{<:AbstractMatrix{<:Number}}) = map(m₀->decompose(m, m₀), ms)

"""
    isintratriangle(
        p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
        vertexes::NTuple{3, Bool}=(true, true, true), edges::NTuple{3, Bool}=(true, true, true), atol::Real=atol, rtol::Real=rtol
    ) -> Bool

Judge whether a point belongs to the interior of a triangle whose vertexes are `p₁`, 'p₂' and `p₃` with the give tolerance. `vertexes` and `edges` define whether the interior should contain the vertexes or edges, respectively.
!!! note
    1. The vertexes are in the order (p₁, p₂, p₃) and the edges are in the order (p₁p₂, p₂p₃, p₃p₁).
    2. The edges do not contain the vertexes.
"""
function isintratriangle(
    p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
    vertexes::NTuple{3, Bool}=(true, true, true), edges::NTuple{3, Bool}=(true, true, true), atol::Real=atol, rtol::Real=rtol
)
    @assert length(p)==length(p₁)==length(p₂)==length(p₃) "isintratriangle error: shape mismatch of input point and triangle."
    @assert 1<length(p)<4 "isintratriangle error: unsupported dimension($(length(p))) of input points."
    x = if length(p) == 2
        p, p₁, p₂, p₃ = SVector{2}(p), SVector{2}(p₁), SVector{2}(p₂), SVector{2}(p₃)
        decompose(p-p₁, p₂-p₁, p₃-p₁)
    else
        p, p₁, p₂, p₃ = SVector{3}(p), SVector{3}(p₁), SVector{3}(p₂), SVector{3}(p₃)
        decompose(p-p₁, p₂-p₁, p₃-p₁)
    end
    x₁_approx_0, x₂_approx_0 = isapprox(x[1], 0.0, atol=atol, rtol=rtol), isapprox(x[2], 0.0, atol=atol, rtol=rtol)
    x₁_approx_1, x₂_approx_1 = isapprox(x[1], 1.0, atol=atol, rtol=rtol), isapprox(x[2], 1.0, atol=atol, rtol=rtol)
    x₁₂_approx_1 = isapprox(x[1]+x[2], 1.0, atol=atol, rtol=rtol)
    x₁_approx_0 && x₂_approx_0 && return vertexes[1]
    x₁_approx_1 && x₂_approx_0 && return vertexes[2]
    x₁_approx_0 && x₂_approx_1 && return vertexes[3]
    x₂_approx_0 && 0<x[1]<1 && return edges[1]
    x₁₂_approx_1 && 0<x[1]<1 && 0<x[2]<1 && return edges[2]
    x₁_approx_0 && 0<x[2]<1 && return edges[3]
    return 0<x[1]<1 && 0<x[2]<1 && x[1]+x[2]<1
end

"""
    isonline(
        p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number};
        ends::Tuple{Bool, Bool}=(true, true), atol::Real=atol, rtol::Real=rtol
    ) -> Bool

Judge whether a point is on a line segment whose end points are `p₁` and `p₂` with the given tolerance. `ends` defines whether the line segment should contain its ends.
"""
function isonline(
    p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number};
    ends::Tuple{Bool, Bool}=(true, true), atol::Real=atol, rtol::Real=rtol
)
    @assert length(p)==length(p₁)==length(p₂) "isonline error: shape mismatch of input point and line segment."
    d₁, d₂, d = distance(p, p₁), distance(p, p₂), distance(p₁, p₂)
    isapprox(d₁, 0.0, atol=atol, rtol=rtol) && return ends[1]
    isapprox(d₂, 0.0, atol=atol, rtol=rtol) && return ends[2]
    return isapprox(d₁+d₂, d, atol=atol, rtol=rtol)
end

"""
    isparallel(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}; atol::Real=atol, rtol::Real=rtol) -> Int

Judge whether two vectors are parallel to each other with the given tolerance, `0` for not parallel, `1` for parallel and `-1` for antiparallel.
"""
function isparallel(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}; atol::Real=atol, rtol::Real=rtol)
    norm₁, norm₂ = norm(v₁), norm(v₂)
    if isapprox(norm₁, 0.0, atol=atol, rtol=rtol) || isapprox(norm₂, 0.0, atol=atol, rtol=rtol)
        result = 1
    elseif length(v₁) == length(v₂)
        temp = dot(v₁, v₂) / norm₁ / norm₂
        result = isapprox(temp, 1, atol=atol, rtol=rtol) ? 1 : isapprox(temp, -1, atol=atol, rtol=rtol) ? -1 : 0
    else
        error("isparallel error: shape mismatch of the input vectors.")
    end
    return result
end

"""
    issubordinate(coordinate::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}; atol::Real=atol, rtol::Real=rtol) -> Bool

Judge whether a coordinate belongs to a lattice defined by `vectors` with the given tolerance.
"""
function issubordinate(coordinate::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}; atol::Real=atol, rtol::Real=rtol)
    @assert length(coordinate)∈(1, 2, 3) "issubordinate error: only 1, 2 and 3 dimensional coordinates are supported."
    @assert length(vectors)∈(1, 2, 3) "issubordinate error: the number of input basis vectors must be 1, 2 or 3."
    fapprox = xi->isapprox(round(xi), xi, atol=atol, rtol=rtol)
    if length(vectors) == 1
        result = mapreduce(fapprox, &, decompose(coordinate, first(vectors)))
    elseif length(vectors) == 2
        result = mapreduce(fapprox, &, decompose(coordinate, first(vectors), last(vectors)))
    else
        v₁, v₂, v₃ = vectors
        result = mapreduce(fapprox, &, decompose(coordinate, v₁, v₂, v₃))
    end
    return result
end

"""
    reciprocals(vectors::AbstractVector{AbstractVector{<:Number}}) -> AbstractVector{<:AbstractVector{<:Number}}

Get the reciprocals dual to the input vectors.
"""
function reciprocals(vectors::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(vectors)<4 "reciprocals error: the number of input vectors should not be greater than 3."
    @assert mapreduce(v->length(v)∈(1, 2, 3), &, vectors, init=true) "reciprocals error: all input vectors must be 1, 2 or 3 dimensional."
    dpi = 2convert(promote_type(Float, eltype(eltype(vectors))), pi)
    if length(vectors) == 1
        v = first(vectors)
        return vector(typeof(vectors), (dpi/mapreduce(vᵢ->vᵢ^2, +, v)*v,))
    elseif length(vectors) == 2
        v₁, v₂ = first(vectors), last(vectors)
        @assert length(v₁)==length(v₂)>1 "reciprocals error: mismatched length of input vectors."
        if length(v₁) == 2
            m = dpi* inv(hcat(SVector{2}(v₁), SVector{2}(v₂)))
            rv₁ = vector(eltype(vectors), Tuple(m[1, :]))
            rv₂ = vector(eltype(vectors), Tuple(m[2, :]))
            return vector(typeof(vectors), (rv₁, rv₂))
        else
            v₁, v₂ = SVector{3}(v₁), SVector{3}(v₂)
            m = dpi * inv(hcat(v₁, v₂, cross(v₁, v₂)))
            rv₁ = vector(eltype(vectors), Tuple(m[1, :]))
            rv₂ = vector(eltype(vectors), Tuple(m[2, :]))
            return vector(typeof(vectors), (rv₁, rv₂))
        end
    elseif length(vectors) == 3
        v₁, v₂, v₃ = vectors
        m = dpi * inv(hcat(SVector{3}(v₁), SVector{3}(v₂), SVector{3}(v₃)))
        rv₁ = vector(eltype(vectors), Tuple(m[1, :]))
        rv₂ = vector(eltype(vectors), Tuple(m[2, :]))
        rv₃ = vector(eltype(vectors), Tuple(m[3, :]))
        return vector(typeof(vectors), (rv₁, rv₂, rv₃))
    end
end
@inline vector(::Type{<:AbstractVector}, v::Tuple) = collect(v)
@inline vector(::Type{<:SVector}, v::Tuple) = SVector(v)

"""
    rotate(vector::AbstractVector{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))) -> AbstractVector{<:Number}
    rotate(cluster::AbstractMatrix{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))) -> AbstractMatrix{<:Number}

Get a rotated vector/cluster of the original one by a certain angle around an axis.

The axis is determined by a point it gets through (`nothing` can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.
!!! note
    1. The result is given by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
    2. Only 2 and 3 dimensional vectors can be rotated.
    3. When the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.
"""
@inline function rotate(vector::AbstractVector{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0)))
    return reshape(rotate(reshape(vector, size(vector, 1), 1), angle; axis=axis), axes(vector, 1))
end
function rotate(cluster::AbstractMatrix{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0)))
    @assert size(cluster, 1)∈(2, 3) "rotate error: only 2 and 3 dimensional vectors can be rotated."
    datatype = promote_type(eltype(cluster), typeof(angle), Float)
    center, theta, phi = (isnothing(axis[1]) ? fill!(similar(cluster, axes(cluster, 1)), zero(datatype)) : axis[1]), axis[2][1], axis[2][2]
    @assert axes(center, 1)==axes(cluster, 1) "rotate error: mismatched shape of the input cluster and the point on axis."
    if length(center) == 2
        @assert isapprox(theta, 0, atol=atol, rtol=rtol) && isapprox(phi, 0, atol=atol, rtol=rtol) "rotate error: both the polar and azimuth of the axis for 2d vectors must be 0."
    end
    cosθ, sinθ = cos(angle), sin(angle)
    k, w = SVector{3, datatype}(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)), zero(MVector{3, datatype})
    result = similar(cluster, datatype, axes(cluster))
    for i in axes(cluster, 2)
        for (l, j) in enumerate(axes(cluster, 1))
            w[l] = cluster[j, i] - center[j]
        end
        inner = dot(k, w)
        outer = cross(k, w)
        for (l, j) in enumerate(axes(cluster, 1))
            result[j, i] = w[l]*cosθ + outer[l]*sinθ + k[l]*inner*(1-cosθ) + center[j]
        end
    end
    return result
end

"""
    translate(cluster::AbstractVector{<:Number}, vector::AbstractVector{<:Number}) -> AbstractVector{<:Number}
    translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number}) -> AbstractMatrix{<:Number}

Get the translated cluster of the original one by a vector.
"""
@inline translate(cluster::AbstractVector{<:Number}, vector::AbstractVector{<:Number}) = cluster + vector 
function translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number})
    @assert axes(cluster, 1)==axes(vector, 1) "translate error: mismatched shape of the input cluster and the translation vector."
    result = similar(cluster, axes(cluster))
    for i in axes(cluster, 2)
        for j in axes(cluster, 1)
            result[j, i] = cluster[j, i] + vector[j]
        end
    end
    return result
end 

"""
    tile(cluster::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations) -> AbstractMatrix{<:Number}
    tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations) -> AbstractMatrix{<:Number}

Tile a supercluster by translations of the input cluster.

Basically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by `vectors` and each set of the translation indices contained in `translations`. When translation vectors are empty, a copy of the original cluster will be returned.
"""
@inline tile(cluster::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations) = tile(reshape(cluster, axes(cluster, 1), 1), vectors, translations)
function tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations)
    length(vectors)==0 && return copy(cluster)
    if length(translations)>0
        @assert length(vectors)==length(first(translations)) "tile error: mismatched shape of input vectors and translations."
        @assert all(v->axes(v, 1)==axes(cluster, 1), vectors) "tile error: mismatched shape of input cluster and vectors."
    end
    datatype = promote_type(eltype(cluster), eltype(eltype(vectors)), Float)
    supercluster = similar(cluster, datatype, axes(cluster, 1), size(cluster, 2)*length(translations))
    disp = similar(cluster, datatype, axes(cluster, 1))
    for (i, translation) in enumerate(translations)
        for k in eachindex(disp)
            disp[k] = zero(datatype)
            for (j, vector) in enumerate(vectors)
                disp[k] += vector[k] * translation[j]
            end
        end
        for (k, j) in enumerate(axes(cluster, 2))
            col = (i-1)*size(cluster, 2) + k
            for row in axes(cluster, 1)
                supercluster[row, col] = cluster[row, j] + disp[row]
            end
        end
    end
    return supercluster
end

"""
    minimumlengths(cluster::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Integer=1; coordination::Integer=12) -> Vector{Float}
    minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Integer=1; coordination::Integer=12) -> Vector{Float}

Use kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.

When the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds across the periodic boundaries, the cluster will be pre-translated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The `coordination` parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain `Inf`. This is a sign that you may need a larger `coordination`. Another situation that `Inf` appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.
"""
@inline function minimumlengths(cluster::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Integer=1; coordination::Integer=12)
    return minimumlengths(reshape(cluster, axes(cluster, 1), 1), vectors, nneighbor; coordination=coordination)
end
function minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Integer=1; coordination::Integer=12)
    @assert nneighbor>=0 "minimumlengths error: input nneighbor must be non negative."
    result = fill(Inf, nneighbor)
    if size(cluster, 2) > 0
        cluster = convert(Matrix{Float}, reshape(cluster, size(cluster)))
        vectors = Vector{Float}[convert(Vector{Float}, reshape(vector, size(vector))) for vector in vectors]
        translations = reshape(collect(product(map(v->-nneighbor:nneighbor, vectors)...)), :)
        for translation in translations
            if length(translation)>0 && mapreduce(≠(0), |, translation)
                remove = map(-, translation)
                filter!(≠(remove), translations)
            end
        end
        supercluster = tile(cluster, vectors, translations)
        for len in flatten(knn(KDTree(supercluster), cluster, nneighbor>0 ? min(nneighbor*coordination, size(supercluster, 2)) : 1, true)[2])
            for (i, minlen) in enumerate(result)
                if isapprox(len, minlen, atol=atol, rtol=rtol)
                    break
                elseif 0.0 < len < minlen
                    nneighbor>0 && (result[i+1:nneighbor] = result[i:nneighbor-1])
                    result[i] = len
                    break
                end
            end
        end
        if nneighbor>0 && mapreduce(isequal(Inf), |, result)
            @warn "minimumlengths warning: larger(>$coordination) coordination or smaller(<$nneighbor) nneighbor may be needed."
        end
    end
    return insert!(result, 1, 0)
end

"""
    Neighbors{K, V<:Number} <: CompositeDict{K, V}

Neighbor vs. bond length maps.
"""
struct Neighbors{K, V<:Number} <: CompositeDict{K, V}
    contents::Dict{K, V}
end
@inline Neighbors(pairs...) = Neighbors(Dict(pairs...))
@inline Neighbors(lengths::AbstractVector{<:Number}, ordinals=0:length(lengths)-1) = Neighbors(ordinal=>length for (ordinal, length) in zip(ordinals, lengths))
@inline Base.max(neighbors::Neighbors) = max(values(neighbors)...)
@inline nneighbor(neighbors::Neighbors{<:Integer}) = max(keys(neighbors)...)
@inline nneighbor(neighbors::Neighbors) = length(neighbors)

"""
    interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Neighbors) -> Vector{Tuple{Int, Int, Int}}

Use kdtree to get the intercluster nearest neighbors.
"""
function interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Neighbors)
    @assert axes(cluster₁, 1)==axes(cluster₂, 1) "interlinks error: mismatched space dimension of input clusters."
    result = Tuple{Int, Int, Int}[]
    length(neighbors)==0 && return result
    indexes₁, indexes₂ = collect(axes(cluster₁, 2)), collect(axes(cluster₂, 2))
    cluster₁, cluster₂ = convert(Matrix{Float}, reshape(cluster₁, size(cluster₁))), convert(Matrix{Float}, reshape(cluster₂, size(cluster₂)))
    exchange = false
    if size(cluster₁, 2) > size(cluster₂, 2)
        cluster₁, cluster₂ = cluster₂, cluster₁
        exchange = true
    end
    for (i, indices) in enumerate(inrange(KDTree(cluster₂), cluster₁, max(neighbors)+atol, true))
        for j in indices
            dist = zero(promote_type(eltype(cluster₁), eltype(cluster₂)))
            for k in axes(cluster₁, 1)
                dist = dist + (cluster₂[k, j]-cluster₁[k, i])^2
            end
            dist = sqrt(dist)
            for (nb, len) in neighbors
                if isapprox(len, dist, atol=atol, rtol=rtol)
                    first, second = exchange ? (j, i) : (i, j)
                    push!(result, (nb, indexes₁[first], indexes₂[second]))
                    break
                end
            end
        end
    end
    return result
end

"""
    Point{N, D<:Number}

A point in a unitcell-described lattice.
"""
struct Point{N, D<:Number}
    site::Int
    rcoordinate::SVector{N, D}
    icoordinate::SVector{N, D}
end
@inline Base.:(==)(point₁::Point, point₂::Point) = ==(efficientoperations, point₁, point₂)
@inline Base.isequal(point₁::Point, point₂::Point) = isequal(efficientoperations, point₁, point₂)
@inline Base.show(io::IO, p::Point) = @printf io "Point(%s, %s, %s)" p.site p.rcoordinate p.icoordinate

"""
    Point(site::Integer, rcoordinate::SVector{N, D}, icoordinate::SVector{N, D}) where {N, D<:Number}
    Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, Val(N))) where N
    Point(site::Integer, rcoordinate::AbstractVector{<:Number}, icoordinate::AbstractVector{<:Number}=zero(SVector{length(rcoordinate), Int}))

Construct a labeled point.
"""
@inline function Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, Val(N))) where N
    datatype = promote_type(eltype(rcoordinate), eltype(icoordinate))
    return Point(site, convert(SVector{N, datatype}, rcoordinate), convert(SVector{N, datatype}, icoordinate))
end
@inline function Point(site::Integer, rcoordinate::AbstractVector{<:Number}, icoordinate::AbstractVector{<:Number}=zero(SVector{length(rcoordinate), eltype(rcoordinate)}))
    datatype = promote_type(eltype(rcoordinate), eltype(icoordinate))
    return Point(site, convert(SVector{length(rcoordinate), datatype}, rcoordinate), convert(SVector{length(icoordinate), datatype}, icoordinate))
end

"""
    dimension(point::Point) -> Int
    dimension(::Type{<:Point{N}}) where N -> Int

Get the spatial dimension of a point.
"""
@inline dimension(point::Point) = dimension(typeof(point))
@inline dimension(::Type{<:Point{N}}) where N = N

"""
    scalartype(::Type{<:Point{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a point.
"""
@inline scalartype(::Type{<:Point{N, D} where N}) where {D<:Number} = D

"""
    isintracell(point::Point) -> Bool

Judge whether a point is intra the unitcell.
"""
@inline isintracell(point::Point) = isapprox(norm(point.icoordinate), 0.0, atol=atol, rtol=rtol)

"""
    Bond{K, P<:Point} <: AbstractVector{P}

A generic bond, which could contains several points.
"""
struct Bond{K, P<:Point} <: AbstractVector{P}
    kind::K
    points::Vector{P}
end
@inline Base.size(bond::Bond) = (length(bond.points),)
@inline Base.firstindex(bond::Bond) = 1
@inline Base.lastindex(bond::Bond) = length(bond.points)
@inline Base.show(io::IO, bond::Bond) = show(io, MIME"text/plain"(), bond)
@inline Base.show(io::IO, ::MIME"text/plain", bond::Bond) = @printf io "Bond(%s, %s)" repr(bond.kind) join(map(string, bond.points), ", ")

"""
    Bond(point::Point)
    Bond(kind, point₁::Point, point₂::Point, points::Point...)

Construct a bond.
"""
@inline Bond(point::Point) = Bond(0, [point])
@inline Bond(kind, point₁::Point, point₂::Point, points::Point...) = Bond(kind, [point₁, point₂, points...])

"""
    getindex(bond::Bond, i::Integer) -> Point

Get the ith point contained in a generic bond.
"""
@inline Base.getindex(bond::Bond, i::Integer) = bond.points[i]

"""
    dimension(bond::Bond) -> Int
    dimension(::Type{<:Bond{K, P} where K}) where {P<:Point} -> Int

Get the space dimension of a concrete bond.
"""
@inline dimension(bond::Bond) = dimension(typeof(bond))
@inline dimension(::Type{<:Bond{K, P} where K}) where {P<:Point} = dimension(P)

"""
    reverse(bond::Bond) -> Bond

Get the reversed bond.
"""
@inline Base.reverse(bond::Bond) = Bond(bond.kind, reverse(bond.points))

"""
    rcoordinate(bond::Bond) -> SVector

Get the rcoordinate of the bond.
"""
@inline function rcoordinate(bond::Bond)
    @assert length(bond)∈(1, 2) "rcoordinate error: not supported for a generic bond which contains $(length(bond)) points."
    length(bond)==1 && return bond[1].rcoordinate
    return bond[2].rcoordinate-bond[1].rcoordinate
end

"""
    icoordinate(bond::Bond) -> SVector

Get the icoordinate of the bond.
"""
@inline function icoordinate(bond::Bond)
    @assert length(bond)∈(1, 2) "icoordinate error: not supported for a generic bond which contains $(length(bond)) points."
    length(bond)==1 && return bond[1].icoordinate
    return bond[2].icoordinate-bond[1].icoordinate
end

"""
    isintracell(bond::Bond) -> Bool

Judge whether a bond is intra the unit cell of a lattice.
"""
@inline isintracell(bond::Bond) = all(point->isapprox(norm(point.icoordinate), 0.0, atol=atol, rtol=rtol), bond)

"""
    AbstractLattice{N, D<:Number, M}

Abstract type of a unitcell-described lattice.

It should have the following contents:
- `name::Symbol`: the name of the lattice
- `coordinates::Matrix{D}`: the coordinates of the lattice
- `vectors::SVector{M, SVector{N, D}}`: the translation vectors of the lattice
"""
abstract type AbstractLattice{N, D<:Number, M} end
@inline contentnames(::Type{<:AbstractLattice}) = (:name, :coordinates, :vectors)
@inline Base.:(==)(lattice₁::AbstractLattice, lattice₂::AbstractLattice) = ==(efficientoperations, lattice₁, lattice₂)
@inline Base.isequal(lattice₁::AbstractLattice, lattice₂::AbstractLattice) = isequal(efficientoperations, lattice₁, lattice₂)
@inline Base.eltype(lattice::AbstractLattice) = eltype(typeof(lattice))
@inline Base.eltype(::Type{<:AbstractLattice{N, D}}) where {N, D<:Number} = SVector{N, D}
@inline Base.iterate(lattice::AbstractLattice, state=1) = state>length(lattice) ? nothing : (lattice[state], state+1)
function Base.show(io::IO, lattice::AbstractLattice)
    @printf io "%s(%s)\n" lattice|>typeof|>nameof getcontent(lattice, :name)
    len = length(lattice)
    if len > 0
        @printf io "  with %s %s:\n" len len==1 ? "point" : "points"
        for i = 1:len
            @printf io "    %s\n" lattice[i]
        end
    end
    len = length(getcontent(lattice, :vectors))
    if len > 0
        @printf io "  with %s translation %s:\n" len len==1 ? "vector" : "vectors"
        for i = 1:len
            @printf io "    %s\n" getcontent(lattice, :vectors)[i]
        end
    end
end

"""
    dimension(lattice::AbstractLattice) -> Int
    dimension(::Type{<:AbstractLattice{N}}) where N -> Int

Get the space dimension of the lattice.
"""
@inline dimension(lattice::AbstractLattice) = dimension(typeof(lattice))
@inline dimension(::Type{<:AbstractLattice{N}}) where N = N

"""
    scalartype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a lattice.
"""
@inline scalartype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number} = D

"""
    length(lattice::AbstractLattice) -> Int

Get the number of points contained in a lattice.
"""
@inline Base.length(lattice::AbstractLattice) = size(getcontent(lattice, :coordinates))[2]

"""
    getindex(lattice::AbstractLattice, i::Integer) -> SVector

Get the ith coordinate.
"""
@inline Base.getindex(lattice::AbstractLattice, i::Integer) = SVector{dimension(lattice), scalartype(lattice)}(ntuple(j->lattice.coordinates[j, i], Val(dimension(lattice))))

"""
    reciprocals(lattice::AbstractLattice) -> Vector{<:SVector}

Get the reciprocal translation vectors of the dual lattice.
"""
@inline reciprocals(lattice::AbstractLattice) = reciprocals(getcontent(lattice, :vectors))

"""
    Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12)

Get the neighbor vs. bond length map of a lattice up to the `nneighbor`th order.
"""
@inline Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12) = Neighbors(minimumlengths(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), nneighbor; coordination=coordination))

"""
    bonds(lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12) -> Vector{Bond{Int, Point{dimension(lattice), scalartype(lattice)}}}
    bonds(lattice::AbstractLattice, neighbors::Neighbors) -> Vector{Bond{keytype(neighbors), Point{dimension(lattice), scalartype(lattice)}}}

Get the required bonds of a lattice.
"""
@inline bonds(lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12) = bonds!(Bond{Int, Point{dimension(lattice), scalartype(lattice)}}[], lattice, nneighbor; coordination=coordination)
@inline bonds(lattice::AbstractLattice, neighbors::Neighbors) = bonds!(Bond{keytype(neighbors), Point{dimension(lattice), scalartype(lattice)}}[], lattice, neighbors)

"""
    bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12)
    bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors) -> typeof(bonds)

Get the required bonds of a lattice and append them to the input bonds.
"""
@inline bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Integer; coordination::Integer=12) = bonds!(bonds, lattice, Neighbors(lattice, nneighbor; coordination=coordination))
function bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors)
    origin = zero(SVector{dimension(lattice), scalartype(lattice)})
    reverse = Dict(length=>order for (order, length) in neighbors)
    haskey(reverse, zero(valtype(neighbors))) && for (i, coordinate) in enumerate(lattice)
        push!(bonds, Bond(reverse[zero(valtype(neighbors))], [Point(i, coordinate, origin)]))
    end
    for (k, index₁, index₂) in interlinks(getcontent(lattice, :coordinates), getcontent(lattice, :coordinates), neighbors)
        if index₂ < index₁
            point₁ = Point(index₁, lattice[index₁], origin)
            point₂ = Point(index₂, lattice[index₂], origin)
            push!(bonds, Bond(k, point₁, point₂))
        end
    end
    nnb = nneighbor(neighbors)
    translations = reshape(collect(product(map(v->-nnb:nnb, getcontent(lattice, :vectors))...)), :)
    for translation in translations
        remove = map(-, translation)
        filter!(≠(remove), translations)
    end
    if length(translations) > 0
        superrcoordinates = tile(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), translations)
        supericoordinates = tile(zero(getcontent(lattice, :coordinates)), getcontent(lattice, :vectors), translations)
        for (k, index₁, index₂) in interlinks(superrcoordinates, getcontent(lattice, :coordinates), neighbors)
            rcoordinate = SVector{dimension(lattice), scalartype(lattice)}(ntuple(j->superrcoordinates[j, index₁], Val(dimension(lattice))))
            icoordinate = SVector{dimension(lattice), scalartype(lattice)}(ntuple(j->supericoordinates[j, index₁], Val(dimension(lattice))))
            point₁ = Point((index₁-1)%length(lattice)+1, rcoordinate, icoordinate)
            point₂ = Point(index₂, lattice[index₂], origin)
            push!(bonds, Bond(k, point₁, point₂))
        end
    end
    return bonds
end

"""
    @recipe plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)

Define the recipe for the visualization of a lattice.
"""
@recipe function plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)
    title --> String(getcontent(lattice, :name))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    framestyle --> :box
    if isa(neighbors, Int) || 0∈keys(neighbors)
        @series begin
            seriestype := :scatter
            markerstrokewidth --> 0
            coordinates = NTuple{dimension(lattice), scalartype(lattice)}[]
            sites = String[]
            for i = 1:length(lattice)
                bond = Bond(Point(i, lattice[i], zero(lattice[i])))
                if filter(bond)
                    push!(coordinates, Tuple(lattice[i]))
                    push!(sites, siteon ? string(i) : "")
                end
            end
            series_annotations := sites
            coordinates
        end
    end
    data = [Float64[] for _ in 1:dimension(lattice)]
    for bond in bonds(lattice, neighbors)
        if length(bond)>0 && filter(bond)
            for i = 1:dimension(lattice)
                for point in bond
                    push!(data[i], point.rcoordinate[i])
                end
                push!(data[i], NaN)
            end
        end
    end
    Tuple(data)
end

"""
    Lattice{N, D<:Number, M} <: AbstractLattice{N, D, M}

Simplest lattice.

A simplest lattice can be constructed from its coordinates and translation vectors.
"""
struct Lattice{N, D<:Number, M} <: AbstractLattice{N, D, M}
    name::Symbol
    coordinates::Matrix{D}
    vectors::SVector{M, SVector{N, D}}
    function Lattice(name::Symbol, coordinates::AbstractMatrix{<:Number}, vectors::SVector{M, <:SVector{N, <:Number}}) where {N, M}
        @assert N==size(coordinates, 1) "Lattice error: shape mismatched."
        datatype = promote_type(Float, eltype(coordinates), eltype(eltype(vectors)))
        coordinates = convert(Matrix{datatype}, coordinates)
        vectors = convert(SVector{M, SVector{N, datatype}}, vectors)
        new{N, datatype, M}(name, coordinates, vectors)
    end
end

"""
    Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    Lattice(coordinates::AbstractVector{<:Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing)

Construct a lattice.
"""
function Lattice(coordinate::NTuple{N, Number}, coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    vectors = isnothing(vectors) ? SVector{0, SVector{N, eltype(eltype(coordinates))}}() : vectorconvert(vectors)
    coordinates = (coordinate, coordinates...)
    coordinates = [coordinates[j][i] for i=1:N, j=1:length(coordinates)]
    return Lattice(name, coordinates, vectors)
end
function Lattice(coordinate::AbstractVector{<:Number}, coordinates::AbstractVector{<:Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing)
    coordinates = hcat(coordinate, coordinates...)
    vectors = isnothing(vectors) ? SVector{0, SVector{size(coordinates)[1], eltype(coordinates)}}() : vectorconvert(vectors)
    return Lattice(name, coordinates, vectors)
end
@inline vectorconvert(vectors::SVector{N, <:SVector}) where N = vectors
@inline vectorconvert(vectors::AbstractVector{<:SVector}) = convert(SVector{length(vectors), SVector{length(eltype(vectors)), eltype(eltype(vectors))}}, vectors)
@inline vectorconvert(vectors::AbstractVector{<:AbstractVector}) = convert(SVector{length(vectors), SVector{length(first(vectors)), eltype(eltype(vectors))}}, vectors)

"""
    Lattice(lattice::AbstractLattice, ranges::OneAtLeast{Int}, boundaries::OneAtLeast{Union{Char, String, Symbol}}=ntuple(i->'O', Val(fieldcount(typeof(ranges)))); mode::Symbol=:nonnegative)
    Lattice(lattice::AbstractLattice, ranges::OneAtLeast{UnitRange{Int}}, boundaries::OneAtLeast{Union{Char, String, Symbol}}=ntuple(i->'O', Val(fieldcount(typeof(ranges)))))

Construct a lattice from the translations of another.
"""
function Lattice(lattice::AbstractLattice, ranges::OneAtLeast{Int}, boundaries::OneAtLeast{Union{Char, String, Symbol}}=ntuple(i->'O', Val(fieldcount(typeof(ranges)))); mode::Symbol=:nonnegative)
    @assert mode∈(:center, :nonnegative) "Lattice error: wrong mode($(repr(mode)))."
    @assert length(ranges)==length(boundaries) "Lattice error: mismatched ranges and boundaries."
    return Lattice(lattice, mode==:center ? map(i->-floor(Int, (i-1)/2):-floor(Int, (i-1)/2)+i-1, ranges) : map(i->0:i-1, ranges), boundaries)
end
function Lattice(lattice::AbstractLattice, ranges::OneAtLeast{UnitRange{Int}}, boundaries::OneAtLeast{Union{Char, String, Symbol}}=ntuple(i->'O', Val(fieldcount(typeof(ranges)))))
    @assert all(map(in(('P', 'O', 'p', 'o', "P", "O", "p", "o", :periodic, :open)), boundaries)) "Lattice error: boundary conditions must be either 'P'/'p'/\"P\"/\"p\"/:periodic for periodic or 'O'/'o'/\"O\"/\"o\"/:open for open."
    @assert length(ranges)==length(boundaries) "Lattice error: mismatched ranges and boundaries."
    boundaries = map(boundary, boundaries)
    name = Symbol(@sprintf "%s%s" getcontent(lattice, :name) join([@sprintf("%s%s%s", boundary=='P' ? "[" : "(", range, boundary=='P' ? "]" : ")") for (range, boundary) in zip(ranges, boundaries)]))
    coordinates = tile(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), product(ranges...))
    vectors = SVector{dimension(lattice), scalartype(lattice)}[]
    for (i, vector) in enumerate(getcontent(lattice, :vectors))
        boundaries[i]=='P' && push!(vectors, vector*length(ranges[i]))
    end
    return Lattice(name, coordinates, vectorconvert(vectors))
end
@inline boundary(b::Union{Char, String}) = first(uppercase(b))
@inline boundary(b::Symbol) = b==:open ? 'O' : 'P'

"""
    expand(momentum::𝕂, reciprocals::AbstractVector{<:AbstractVector{<:Number}}) -> eltype(reciprocals)

Expand the momentum from integral values to real values with the given reciprocals.
"""
@inline function expand(momentum::𝕂, reciprocals::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(momentum)==length(reciprocals) "expand error: mismatched momentum and reciprocals."
    vs, ps = values(momentum), periods(momentum)
    result = zero(first(reciprocals))
    for (i, index) in enumerate(eachindex(reciprocals))
        result += reciprocals[index] * (vs[i]//ps[i])
    end
    return result
end

"""
    𝕂¹{N}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where N

Construct a 1d quantum momentum by the coordinates.
"""
function 𝕂¹{N}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where N
    @assert length(reciprocals)==1 "𝕂¹ error: mismatched length of reciprocals."
    k = decompose(momentum, first(reciprocals)) * N
    i = round(Int, k)
    @assert isapprox(i, k; atol=atol, rtol=rtol) "𝕂¹ error: input momentum not on grid."
    return 𝕂¹{N}(i)
end

"""
    𝕂²{N₁, N₂}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where {N₁, N₂}

Construct a 2d quantum momentum by the coordinates.
"""
function 𝕂²{N₁, N₂}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where {N₁, N₂}
    @assert length(reciprocals)==2 "𝕂² error: mismatched length of reciprocals."
    k₁, k₂ = decompose(momentum, first(reciprocals), last(reciprocals))
    k₁, k₂ = k₁*N₁, k₂*N₂
    i₁, i₂ = round(Int, k₁), round(Int, k₂)
    @assert isapprox(i₁, k₁; atol=atol, rtol=rtol) && isapprox(i₂, k₂; atol=atol, rtol=rtol) "𝕂² error: input momentum not on grid."
    return 𝕂²{N₁, N₂}(i₁, i₂)
end

"""
    𝕂³{N₁, N₂, N₃}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where {N₁, N₂, N₃}

Construct a 3d quantum momentum by the coordinates.
"""
function 𝕂³{N₁, N₂, N₃}(momentum::AbstractVector{<:Number}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}; atol=atol, rtol=rtol) where {N₁, N₂, N₃}
    @assert length(reciprocals)==3 "𝕂³ error: mismatched length of reciprocals."
    v₁, v₂, v₃ = reciprocals
    k₁, k₂, k₃ = decompose(momentum, v₁, v₂, v₃)
    k₁, k₂, k₃ = k₁*N₁, k₂*N₂, k₃*N₃
    i₁, i₂, i₃ = round(Int, k₁), round(Int, k₂), round(Int, k₃)
    @assert isapprox(i₁, k₁; atol=atol, rtol=rtol) && isapprox(i₂, k₂; atol=atol, rtol=rtol) && isapprox(i₃, k₃; atol=atol, rtol=rtol) "𝕂³ error: input momentum not on grid."
    return 𝕂³{N₁, N₂, N₃}(i₁, i₂, i₃)
end

"""
    ReciprocalSpace{K, P<:SVector} <: VectorSpace{P}

Abstract type of reciprocal spaces.
"""
abstract type ReciprocalSpace{K, P<:SVector} <: VectorSpace{P} end

"""
    scalartype(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector}

Get the data type of the coordinates of a reciprocal space.
"""
@inline scalartype(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector} = eltype(P)

"""
    dimension(reciprocalspace::ReciprocalSpace) -> Int
    dimension(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector} -> Int

Get the spatial dimension of a reciprocal space.
"""
@inline dimension(reciprocalspace::ReciprocalSpace) = dimension(typeof(reciprocalspace))
@inline dimension(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector} = length(P)

"""
    label(reciprocalspace::ReciprocalSpace) -> Symbol
    label(::Type{<:ReciprocalSpace{K}}) where K -> Symbol

Get the label of a reciprocal space.
"""
@inline label(reciprocalspace::ReciprocalSpace) = label(typeof(reciprocalspace))
@inline label(::Type{<:ReciprocalSpace{K}}) where K = K

"""
    label(reciprocalspace::ReciprocalSpace, i::Integer) -> String
    label(::Type{<:ReciprocalSpace{K}}, i::Integer) where K -> String

Get the label of the ith axis of a reciprocal space.
"""
@inline label(reciprocalspace::ReciprocalSpace, i::Integer) = label(typeof(reciprocalspace), i)
@inline label(::Type{<:ReciprocalSpace{K}}, i::Integer) where K = string(K, subscript(i))

"""
    FractionalReciprocalSpace{K, P<:SVector, N} <: ReciprocalSpace{K, P}

Abstract type of reciprocal spaces with fractional coordinates.
"""
abstract type FractionalReciprocalSpace{K, N, P<:SVector} <: ReciprocalSpace{K, P} end

"""
    reciprocals(reciprocalspace::FractionalReciprocalSpace) -> AbstractVector{<:AbstractVector{<:Number}}

Get the reciprocal lattice vectors of a reciprocal space with fractional coordinates.
"""
@inline reciprocals(reciprocalspace::FractionalReciprocalSpace) = getcontent(reciprocalspace, :reciprocals)

"""
    @recipe plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true)

Define the recipe for the visualization of a reciprocal space with fractional coordinates.

When `fractional` is `true`, the fractional coordinates will be plotted. Otherwise the Cartesian coordinates will be plotted.
"""
@recipe function plot(reciprocalspace::FractionalReciprocalSpace{K, N}; fractional=false, autolims=true) where {K, N}
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    aspect_ratio --> :equal
    legend --> false
    minorgrid --> true
    framestyle --> :box
    if fractional
        coordinates = map(NTuple{N, scalartype(reciprocalspace)}, fractionals(reciprocalspace))
        N>0 && (xlabel-->label(reciprocalspace, 1); autolims || (xlims-->extrema(v->v[1], coordinates)))
        N>1 && (ylabel-->label(reciprocalspace, 2); autolims || (ylims-->extrema(v->v[2], coordinates)))
        N>2 && (zlabel-->label(reciprocalspace, 3); autolims || (zlims-->extrema(v->v[3], coordinates)))
        coordinates
    else
        dim = dimension(reciprocalspace)
        coordinates = map(NTuple{dim, scalartype(reciprocalspace)}, reciprocalspace)
        dim>0 && (xlabel-->label(reciprocalspace, 1); autolims || (xlims-->extrema(v->v[1], coordinates)))
        dim>1 && (ylabel-->label(reciprocalspace, 2); autolims || (ylims-->extrema(v->v[2], coordinates)))
        dim>2 && (zlabel-->label(reciprocalspace, 3); autolims || (zlims-->extrema(v->v[3], coordinates)))
        coordinates
    end
end

"""
    ProductedReciprocalSpace{K, N, P<:SVector} <: FractionalReciprocalSpace{K, N, P}

Abstract type of reciprocal spaces whose fractional coordinates are composed of the direct product of several ranges.
"""
abstract type ProductedReciprocalSpace{K, N, P<:SVector} <: FractionalReciprocalSpace{K, N, P} end
@inline VectorSpaceStyle(::Type{<:ProductedReciprocalSpace}) = VectorSpaceDirectProducted(:backward)

"""
    fractionals(reciprocalspace::ProductedReciprocalSpace{K, N}) where {K, N} -> AbstractVector{SVector{N, scalartype(reciprocalspace)}}

Get the fractional coordinates of a reciprocal space.
"""
@inline function fractionals(reciprocalspace::ProductedReciprocalSpace{K, N}) where {K, N}
    return DirectProductedVectorSpace{:backward, SVector{N, scalartype(reciprocalspace)}}(ntuple(i->range(reciprocalspace, i), Val(N)))
end

"""
    BrillouinZone{K, P<:𝕂, N, S<:SVector} <: ProductedReciprocalSpace{K, N, S}

Brillouin zone of a lattice.
"""
struct BrillouinZone{K, P<:𝕂, N, S<:SVector} <: ProductedReciprocalSpace{K, N, S}
    reciprocals::SVector{N, S}
    function BrillouinZone{K, P}(reciprocals::SVector{N, <:SVector}) where {K, P<:𝕂, N}
        @assert isa(K, Symbol) "BrillouinZone error: K must be a Symbol."
        @assert rank(P)==N "BrillouinZone error: mismatched momentum and reciprocals."
        new{K, P, N, eltype(reciprocals)}(reciprocals)
    end
end
@inline shape(::BrillouinZone{K, P}) where {K, P<:𝕂} = shape(Momenta(P))
@inline Base.convert(::Type{<:SVector}, index::CartesianIndex, brillouinzone::BrillouinZone{K, P}) where {K, P<:𝕂} = expand(Momenta(P)[index], brillouinzone.reciprocals)

"""
    BrillouinZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, nk)
    BrillouinZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, nk) where K
    BrillouinZone(::Type{P}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}) where {P<:𝕂}
    BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}) where {K, P<:𝕂}

Construct a Brillouin zone.
"""
@inline BrillouinZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, nk) = BrillouinZone{:k}(reciprocals, nk)
@inline BrillouinZone(::Type{P}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}) where {P<:𝕂} = BrillouinZone{:k}(P, reciprocals)
@inline function BrillouinZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, nk) where K
    isa(nk, Integer) && (nk = map(v->nk, reciprocals))
    @assert length(nk)==length(reciprocals) "BrillouinZone error: mismatched number of reciprocals and periods of momenta."
    length(reciprocals)==1 && return BrillouinZone{K}(𝕂¹{nk[1]}, reciprocals)
    length(reciprocals)==2 && return BrillouinZone{K}(𝕂²{nk[1], nk[2]}, reciprocals)
    length(reciprocals)==3 && return BrillouinZone{K}(𝕂³{nk[1], nk[2], nk[3]}, reciprocals)
    error("BrillouinZone error: only 1d, 2d and 3d are supported.")
end
@inline BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector{<:AbstractVector{<:Number}}) where {K, P<:𝕂} = BrillouinZone{K, P}(vectorconvert(reciprocals))

"""
    hash(brillouinzone::BrillouinZone, h::UInt)

Hash a Brillouin zone.
"""
@inline Base.hash(brillouinzone::BrillouinZone, h::UInt) = hash((Tuple(brillouinzone.reciprocals), periods(keytype(brillouinzone))), h)

"""
    keytype(brillouinzone::BrillouinZone)
    keytype(::Type{<:BrillouinZone{K, P} where K}) where {P<:𝕂}

Get the keytype of a Brillouin zone, which is defined to be the type of the Abelian quantum number representing a momentum in the Brillouin zone.
"""
@inline Base.keytype(brillouinzone::BrillouinZone) = keytype(typeof(brillouinzone))
@inline Base.keytype(::Type{<:BrillouinZone{K, P} where K}) where {P<:𝕂} = P

"""
    keys(::BrillouinZone{K, P}) where {K, P<:𝕂} -> Momenta{P}

Get the keys of a Brillouin zone, which is defined to be the complete allowed set of Abelian quantum numbers representing the momenta in the Brillouin zone.
"""
@inline Base.keys(::BrillouinZone{K, P}) where {K, P<:𝕂} = Momenta(P)

"""
    getindex(brillouinzone::BrillouinZone{K, P}, index::P) where {K, P<:𝕂} -> eltype(brillouinzone)

Get the coordinates of a momentum in a Brillouin zone by its corresponding Abelian quantum number.
"""
@inline Base.getindex(brillouinzone::BrillouinZone{K, P}, index::P) where {K, P<:𝕂} = expand(index, brillouinzone.reciprocals)

"""
    range(brillouinzone::BrillouinZone, i::Integer) -> StepRangeLen{Float64}

Get the ith axis range of the fractional coordinates of a Brillouin zone.
"""
@inline function Base.range(brillouinzone::BrillouinZone, i::Integer)
    n = period(keytype(brillouinzone), i)
    return range(0, (n-1)/n, n)
end

"""
    volume(brillouinzone::BrillouinZone) -> Number

Get the volume of a Brillouin zone.
"""
@inline volume(brillouinzone::BrillouinZone) = volume(brillouinzone.reciprocals)

"""
    ReciprocalZone{K, N, S<:SVector, V<:Number} <: ProductedReciprocalSpace{K, N, S}

A zone in the reciprocal space.
"""
struct ReciprocalZone{K, N, S<:SVector, V<:Number} <: ProductedReciprocalSpace{K, N, S}
    reciprocals::SVector{N, S}
    bounds::NTuple{N, Segment{V}}
    volume::V
    function ReciprocalZone{K}(reciprocals::SVector{N, <:SVector}, bounds::NTuple{N, Segment}) where {K, N}
        @assert isa(K, Symbol) "ReciprocalZone error: K must be a Symbol."
        ratio = one(promote_type(eltype(eltype(reciprocals)), eltype(eltype(bounds))))
        for i = 1:length(bounds)
            ratio = ratio * (bounds[i].stop-bounds[i].start)
        end
        v = ratio * volume(reciprocals)
        new{K, N, eltype(reciprocals), typeof(v)}(reciprocals, bounds, v)
    end
end
@inline shape(reciprocalzone::ReciprocalZone) = map(bound->1:length(bound), reciprocalzone.bounds)
@inline function Base.convert(::Type{<:SVector}, index::CartesianIndex{N}, reciprocalzone::ReciprocalZone{K, N}) where {K, N}
    result = zero(eltype(reciprocalzone))
    for (reciprocal, bound, i) in zip(reciprocalzone.reciprocals, reciprocalzone.bounds, index.I)
        result += reciprocal * bound[i]
    end
    return result
end

"""
    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}; length=100, ends=(true, false))
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}; length=100, ends=(true, false)) where K

    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false))
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false)) where K

    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::AbstractVector{<:Pair{<:Number, <:Number}}; length=100, ends=(true, false))
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::AbstractVector{<:Pair{<:Number, <:Number}}; length=100, ends=(true, false)) where K

    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::OneAtLeast{Pair{<:Number, <:Number}}; length=100, ends=(true, false))
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::OneAtLeast{Pair{<:Number, <:Number}}; length=100, ends=(true, false)) where K

Construct a rectangular zone in the reciprocal space.
"""
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals); length=length, ends=ends)
end
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bound::Pair{<:Number, <:Number}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals), bound, bounds...; length=length, ends=ends)
end
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::AbstractVector{<:Pair{<:Number, <:Number}}; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals), Tuple(bounds); length=length, ends=ends)
end
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::OneAtLeast{Pair{<:Number, <:Number}}; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals), bounds; length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}; length=100, ends=(true, false)) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), ntuple(i->-1//2=>1//2, Val(Base.length(reciprocals))); length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bound::Pair{<:Number, <:Number}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false)) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), (bound, bounds...); length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::AbstractVector{<:Pair{<:Number, <:Number}}; length=100, ends=(true, false)) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), Tuple(bounds); length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, bounds::OneAtLeast{Pair{<:Number, <:Number}}; length=100, ends=(true, false)) where K
    N = fieldcount(typeof(bounds))
    isa(length, Integer) && (length = ntuple(i->length, Val(N)))
    isa(ends, NTuple{2, Bool}) && (ends = ntuple(i->ends, Val(N)))
    @assert Base.length(reciprocals)==N "ReciprocalZone error: the numbers of the input reciprocals($(Base.length(reciprocals))) and bounds($N) do not match."
    @assert Base.length(length)==N "ReciprocalZone error: the number of length should be $N if it is not an integer."
    @assert Base.length(ends)==N "ReciprocalZone error: the number of ends should be $N if it is not `NTuple{2, Bool}`."
    segments = ntuple(i->Segment(bounds[i].first, bounds[i].second, length[i]; ends=ends[i]), Val(N))
    return ReciprocalZone{K}(vectorconvert(reciprocals), segments)
end

"""
    ReciprocalZone(reciprocalspace::ProductedReciprocalSpace)

Construct a reciprocal zone from a reciprocal space whose fractional coordinates are composed of the direct product of several ranges.
"""
@inline function ReciprocalZone(reciprocalspace::ProductedReciprocalSpace{K, N}) where {K, N}
    return ReciprocalZone{K}(
        reciprocals(reciprocalspace),
        ntuple(i->(bound=range(reciprocalspace, i); Segment(first(bound), last(bound)+step(bound), length(bound); ends=(true, false))), Val(N))
    )
end

"""
    range(reciprocalzone::ReciprocalZone, i::Integer) -> StepRangeLen{Float64}

Get the ith axis range of the fractional coordinates of a reciprocal zone.
"""
@inline function Base.range(reciprocalzone::ReciprocalZone, i::Integer)
    bound = reciprocalzone.bounds[i]
    return range(first(bound), last(bound), length(bound))
end

"""
    shrink(reciprocalzone::ReciprocalZone{K}, ranges::Vararg{OrdinalRange{<:Integer}, N}) where {K, N} -> ReciprocalZone

Shrink a reciprocal zone.
"""
function shrink(reciprocalzone::ReciprocalZone{K}, ranges::Vararg{OrdinalRange{<:Integer}, N}) where {K, N}
    @assert length(ranges)==length(reciprocalzone.reciprocals) "shrink error: mismatched number of ranges and reciprocals."
    return ReciprocalZone{K}(reciprocalzone.reciprocals, ntuple(i->reciprocalzone.bounds[i][ranges[i]], Val(N)))
end

"""
    volume(reciprocalzone::ReciprocalZone) -> Number

Get the volume of a reciprocal zone.
"""
@inline volume(reciprocalzone::ReciprocalZone) = reciprocalzone.volume

"""
    ReciprocalScatter{K, N, S<:SVector, V<:SVector{N}} <: FractionalReciprocalSpace{K, N, S}

A set of scatter points in the reciprocal space.
"""
struct ReciprocalScatter{K, N, S<:SVector, V<:SVector{N}} <: FractionalReciprocalSpace{K, N, S}
    reciprocals::SVector{N, S}
    coordinates::Vector{V}
    function ReciprocalScatter{K}(reciprocals::SVector{N, <:SVector}, coordinates::Vector{<:SVector{N}}) where {K, N}
        @assert isa(K, Symbol) "ReciprocalScatter error: K must be a Symbol."
        new{K, N, eltype(reciprocals), eltype(coordinates)}(reciprocals, coordinates)
    end
end
@inline Base.length(reciprocalscatter::ReciprocalScatter) = length(reciprocalscatter.coordinates)
@inline Base.getindex(reciprocalscatter::ReciprocalScatter, i::Integer) = mapreduce(*, +, reciprocalscatter.reciprocals, reciprocalscatter.coordinates[i])

"""
    ReciprocalScatter(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, coordinates::AbstractVector{<:AbstractVector{<:Number}})
    ReciprocalScatter{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, coordinates::AbstractVector{<:AbstractVector{<:Number}}) where K

Construct a set of reciprocal scatter points.
"""
@inline ReciprocalScatter(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, coordinates::AbstractVector{<:AbstractVector{<:Number}}) = ReciprocalScatter{:k}(reciprocals, coordinates)
function ReciprocalScatter{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, coordinates::AbstractVector{<:AbstractVector{<:Number}}) where K
    reciprocals = vectorconvert(reciprocals)
    coordinates = collect(SVector{length(reciprocals), eltype(eltype(reciprocals))}, coordinates)
    return ReciprocalScatter{K}(reciprocals, coordinates)
end

"""
    ReciprocalScatter(reciprocalspace::FractionalReciprocalSpace)

Construct a set of reciprocal scatter points from a reciprocal space with fractional coordinates.
"""
@inline ReciprocalScatter(reciprocalspace::FractionalReciprocalSpace) = ReciprocalScatter{label(reciprocalspace)}(reciprocals(reciprocalspace), fractionals(reciprocalspace))

"""
    fractionals(reciprocalscatter::ReciprocalScatter) -> Vector{<:SVector}

Get the fractional coordinates of a reciprocal space.
"""
@inline fractionals(reciprocalscatter::ReciprocalScatter) = reciprocalscatter.coordinates

"""
    ReciprocalPath{K, S<:SVector, N, R} <: ReciprocalSpace{K, S}

A path in the reciprocal space.
"""
struct ReciprocalPath{K, S<:SVector, N, R} <: ReciprocalSpace{K, S}
    contents::NTuple{N, Vector{S}}
    labels::NTuple{N, Pair{R, R}}
    function ReciprocalPath{K}(contents::NTuple{N, AbstractVector{<:AbstractVector{<:Number}}}, labels::NTuple{N, Pair{R, R}}) where {K, N, R}
        @assert isa(K, Symbol) "ReciprocalPath error: K must be a Symbol."
        new{K, eltype(eltype(contents)), N, R}(map(vectorconvert, contents), labels)
    end
end
@inline ReciprocalPath(contents::Tuple{Vararg{AbstractVector{<:AbstractVector{<:Number}}}}, labels::Tuple{Vararg{Pair}}) = ReciprocalPath{:k}(contents, labels)
@inline contentnames(::Type{<:ReciprocalPath}) = (:contents, :labels)
@inline VectorSpaceStyle(::Type{<:ReciprocalPath}) = VectorSpaceDirectSummed()

"""
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing)
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing) where K

    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, points::OneOrMore{Number}...; labels=points, length=100, ends=nothing)
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, points::OneOrMore{Number}...; labels=points, length=100, ends=nothing) where K

    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, length=100, ends=nothing)
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, length=100, ends=nothing) where K

Construct a path in the reciprocal space.

!!! note
    1) For connected segments,
       * when `length` is an integer, it specifies the length of each segment except for the last whose length will be `length+1`;
       * when `ends` is `nothing`, the start point will be included while the end point will be not for each segment except for the last both points of which will be included.

    2) For disconnected segments, they can be partitioned into several connected parts, and the rules for connected segments apply for each of such connected parts.

    With the above rules, all the points along the assigned path will be counted once and only once with the largest homogeneity.
"""
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing)
    return ReciprocalPath(reciprocals, contents.points...; labels=contents.labels, length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, point₁::OneOrMore{Number}, point₂::OneOrMore{Number}, points::OneOrMore{Number}...; labels=(point₁, point₂, points...), length=100, ends=nothing)
    return ReciprocalPath(reciprocals, points2segments((point₁, point₂, points...))...; labels=points2segments(Tuple(labels)), length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, length=100, ends=nothing)
    return ReciprocalPath{:k}(reciprocals, segments...; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing) where K
    return ReciprocalPath{K}(reciprocals, contents.points...; labels=contents.labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, point₁::OneOrMore{Number}, point₂::OneOrMore{Number}, points::OneOrMore{Number}...; labels=(point₁, point₂, points...), length=100, ends=nothing) where K
    return ReciprocalPath{K}(reciprocals, points2segments((point₁, point₂, points...))...; labels=points2segments(labels), length=length, ends=ends)
end
@inline @inbounds function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector{<:Number}}, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, length=100, ends=nothing) where K
    segments = map(segment->OneOrMore(segment.first)=>OneOrMore(segment.second), segments)
    N = fieldcount(typeof(first(first(segments))))
    @assert N>0 "ReciprocalPath error: empty path."
    @assert all(segment->Base.length(segment.first)==Base.length(segment.second)==N, segments) "ReciprocalPath error: mismatched dimension of input segments."
    @assert Base.length(reciprocals)==N "ReciprocalPath error: mismatched input reciprocals and segments."
    M = fieldcount(typeof(segments))
    isa(length, Integer) && (length = ntuple(i->i<M ? (segments[i].second==segments[i+1].first ? length : length+1) : length+1, Val(M)))
    isnothing(ends) && (ends = ntuple(i->i<M ? (true, segments[i].second==segments[i+1].first ? false : true) : (true, true), Val(M)))
    @assert Base.length(length)==M "ReciprocalPath error: the number of length should be $M if it is not an integer."
    @assert Base.length(ends)==M "ReciprocalPath error: the number of ends should be $M if it is not `nothing`."
    @assert all((segments[i+1].first==segments[i].second)==(labels[i+1].first==labels[i].second) for i = 1:M-1) "ReciprocalPath error: labels are not consistent with segments."
    reciprocals = vectorconvert(reciprocals)
    segments = map((segment, length, ends)->Segment(segment.first, segment.second, length; ends=ends), segments, length, ends)
    contents = map(segment->zeros(SVector{Base.length(eltype(reciprocals)), promote_type(eltype(eltype(reciprocals)), eltype(eltype(eltype(segments))))}, Base.length(segment)), segments)
    for (m, segment) in enumerate(segments)
        for (n, index) in enumerate(segment)
            for i = 1:Base.length(index)
                contents[m][n] += reciprocals[i]*index[i]
            end
        end
    end
    return ReciprocalPath{K}(contents, homogenization(Tuple(labels)))
end
@inline points2segments(points::NTuple{M, Any}) where M = ntuple(i->points[i]=>points[i+1], Val(M-1))
@inline homogenization(labels::OneAtLeast{Pair{<:Number, <:Number}}) = promote(labels...)
@inline function homogenization(labels::OneAtLeast{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}) where N
    D = promote_type(concatenate(map(fieldtypes, concatenate(map(fieldtypes, fieldtypes(typeof(labels)))...))...)...)
    return ntuple(i->NTuple{N, D}(labels[i].first)=>NTuple{N, D}(labels[i].second), Val(fieldcount(typeof(labels))))
end
@inline homogenization(labels::OneAtLeast{Pair}) = labels

"""
    step(path::ReciprocalPath, i::Integer) -> scalartype(path)

Get the step between the ith and the (i+1)th points in the path.
"""
@inline Base.step(path::ReciprocalPath, i::Integer) = distance(path[i], path[i+1])

"""
    distance(path::ReciprocalPath) -> scalartype(path)
    distance(path::ReciprocalPath, i::Integer) -> scalartype(path)
    distance(path::ReciprocalPath, i::Integer, j::Integer) -> scalartype(path)

Get the distance
1) of the total path,
2) from the start point to the ith point in the path,
3) from the ith point to the jth point in the path (when i is greater than j, the value is negative).
"""
@inline distance(path::ReciprocalPath) = distance(path, 1, length(path))
@inline distance(path::ReciprocalPath, i::Integer) = distance(path, 1, i)
function distance(path::ReciprocalPath, i::Integer, j::Integer)
    i==j && return zero(scalartype(path))
    i>j && return -distance(path, j, i)
    dimsum = cumsum(map(length, path.contents))
    i₁ = searchsortedfirst(dimsum, i)
    j₁ = searchsortedfirst(dimsum, j)
    i₂ = i₁>1 ? (i-dimsum[i₁-1]) : i
    j₂ = j₁>1 ? (j-dimsum[j₁-1]) : j
    i₁==j₁ && return distance(path.contents[i₁][i₂], path.contents[j₁][j₂])
    result = distance(path.contents[i₁][i₂], path.contents[i₁][end]) + distance(path.contents[j₁][1], path.contents[j₁][j₂])
    for l = (i₁+1):(j₁-1)
        path.labels[l-1].second==path.labels[l].first && (result += distance(path.contents[l-1][end], path.contents[l][1]))
        result += distance(path.contents[l][1], path.contents[l][end])
    end
    path.labels[j₁-1].second==path.labels[j₁].first && (result += distance(path.contents[j₁-1][end], path.contents[j₁][1]))
    return result
end

"""
    ticks(path::ReciprocalPath) -> Tuple{Vector{scalartype(path)}, Vector{String}}

Get the position-label pairs of the ticks of a path.
"""
function ticks(path::ReciprocalPath)
    result = (scalartype(path)[], String[])
    d = zero(scalartype(path))
    for i = 1:length(path.contents)
        if i==1
            push!(result[2], string(path.labels[i].first))
        else
            previous, current = path.labels[i-1].second, path.labels[i].first
            previous==current && (d += distance(path.contents[i-1][end], path.contents[i][1]))
            push!(result[2], previous==current ? string(previous) : string(previous, " / ", current))
        end
        push!(result[1], d)
        d += distance(path.contents[i][1], path.contents[i][end])
    end
    push!(result[1], d)
    push!(result[2], string(path.labels[end].second))
    return result
end

"""
    @recipe plot(path::ReciprocalPath)

Define the recipe for the visualization of a reciprocal path.
"""
@recipe function plot(path::ReciprocalPath)
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    legend --> false
    aspect_ratio --> :equal
    minorgrid --> true
    framestyle --> :box
    dim = dimension(path)
    dim>0 && (xlabel --> label(path, 1))
    dim>1 && (ylabel --> label(path, 2))
    dim>2 && (zlabel --> label(path, 3))
    x, y, text = [], [], []
    for i = 1:length(path.contents)
        if length(path.contents[i])>0
            push!(x, path.contents[i][1][1])
            push!(y, path.contents[i][1][2])
            push!(text, string(path.labels[i].first))
            j = i%length(path.labels)+1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                push!(x, path.contents[i][end][1])
                push!(y, path.contents[i][end][2])
                push!(text, string(path.labels[i].second))
            end
        end
    end
    annotation := (x, y, text)
    coordinates = map(Tuple, path)
end

"""
    selectpath(brillouinzone::BrillouinZone, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; ends=nothing, atol::Real=atol, rtol::Real=rtol) -> Tuple(ReciprocalPath, Vector{Int})
    selectpath(brillouinzone::BrillouinZone, points::OneOrMore{Number}...; labels=points, ends=nothing, atol::Real=atol, rtol::Real=rtol) -> Tuple(ReciprocalPath, Vector{Int})
    selectpath(brillouinzone::BrillouinZone, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) -> Tuple(ReciprocalPath, Vector{Int})

Select a path from a `BrillouinZone`. Return a `ReciprocalPath` and the positions of the equivalent points in the `BrillouinZone`.

!!! note
    1) For connected segments, the start point will be included while the stop point will be not for each segment except for the last both points of which will be included if `ends` is `nothing`.
    2) For disconnected segments, they can be partitioned into several connected parts, and the rule for connected segments applies for each of such connected parts.

    With the above rules, all the points along the assigned path will be counted once and only once.
"""
@inline function selectpath(brillouinzone::BrillouinZone, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; ends=nothing, atol::Real=atol, rtol::Real=rtol)
    return selectpath(brillouinzone, contents.points...; labels=contents.labels, ends=ends, atol=atol, rtol=rtol)
end
@inline function selectpath(brillouinzone::BrillouinZone, point₁::OneOrMore{Number}, point₂::OneOrMore{Number}, points::OneOrMore{Number}...; labels=(point₁, point₂, points...), ends=nothing, atol::Real=atol, rtol::Real=rtol)
    return selectpath(brillouinzone, points2segments((point₁, point₂, points...))...; labels=points2segments(Tuple(labels)), ends=ends, atol=atol, rtol=rtol)
end
@inline function selectpath(brillouinzone::BrillouinZone{K}, segments::Pair{<:OneOrMore{Number}, <:OneOrMore{Number}}...; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) where K
    N, M = fieldcount(typeof(first(first(segments)))), fieldcount(typeof(segments))
    @assert length(brillouinzone.reciprocals)==N "selectpath error: mismatched number of reciprocals ($(length(brillouinzone.reciprocals)) v.s. $N)."
    isnothing(ends) && (ends = ntuple(i->i<M ? (true, segments[i].second==segments[i+1].first ? false : true) : (true, true), Val(M)))
    @assert length(ends)==M "selectpath error: the number of ends should be $M if it is not `nothing`"
    contents, indexes = ntuple(i->eltype(brillouinzone)[], Val(M)), Int[]
    for (n, segment) in enumerate(segments)
        start = trunc.(Int, segment.first)
        stop = trunc.(Int, segment.second)
        mesh = StepRange{Int, Int}[]
        for i in 1:N
            if stop[i] >= start[i]
                step = 1
            elseif stop[i] < start[i]
                step = -1
            end
            push!(mesh, start[i]:step:stop[i])
        end
        translations = [mapreduce(*, +, translation, brillouinzone.reciprocals) for translation in product(mesh...)]
        start = mapreduce(*, +, segment.first, brillouinzone.reciprocals)
        stop = mapreduce(*, +, segment.second, brillouinzone.reciprocals)
        positions, distances = Int[], scalartype(brillouinzone)[]
        for (pos, k) in enumerate(brillouinzone)
            for translation in translations
                coord = k + translation
                if isonline(coord, start, stop; ends=ends[n], atol=atol, rtol=rtol)
                    push!(contents[n], coord)
                    push!(positions, pos)
                    push!(distances, distance(coord, start))
                end
            end
        end
        permutation = sortperm(distances)
        permute!(contents[n], permutation)
        positions = positions[permutation]
        append!(indexes, positions)
    end
    return ReciprocalPath{K}(contents, Tuple(labels)), indexes
end

const linemap = Dict(
    "Γ"=>0, "X"=>1//2, "Γ₁"=>0, "Γ₂"=>1, "X₁"=>1//2, "X₂"=>-1//2
)
"""
    line"P₁-P₂-P₃-..."

Construct a tuple of start-stop point pairs for the one dimensional reciprocal space.
"""
macro line_str(str::String)
    points = split(str, "-")
    @assert length(points)>1 "@line_str error: too few points."
    return (points=ntuple(i->linemap[points[i]], length(points)), labels=ntuple(i->points[i], length(points)))
end

const rectanglemap = Dict(
    "Γ"=>(0, 0), "X"=>(1//2, 0), "Y"=>(0, 1//2), "M"=>(1//2, 1//2),
    "X₁"=>(1//2, 0), "X₂"=>(-1//2, 0), "Y₁"=>(0, 1//2), "Y₂"=>(0, -1//2),
    "M₁"=>(1//2, 1//2), "M₂"=>(-1//2, 1//2), "M₃"=>(-1//2, -1//2), "M₄"=>(1//2, -1//2)
)
"""
    rectangle"P₁-P₂-P₃-..."

Construct a tuple of start-stop point pairs for the rectangular reciprocal space.
"""
macro rectangle_str(str::String)
    points = split(replace(str, " "=>""), "-")
    @assert length(points)>1 "@rectangle_str error: too few points."
    return (points=ntuple(i->rectanglemap[points[i]], length(points)), labels=ntuple(i->points[i], length(points)))
end

# High-symmetric point in the hexagonal Brillouin zone when the angle between the reciprocal vectors is 120°
const hexagon120°map = Dict(
    "Γ"=>(0, 0), "K"=>(2//3, 1//3), "M"=>(1//2, 1//2),
    "K₁"=>(2//3, 1//3), "K₂"=>(1//3, 2//3), "K₃"=>(1//3, -1//3), "K₄"=>(-2//3, -1//3), "K₅"=>(-1//3, -2//3), "K₆"=>(-1//3, 1//3),
    "M₁"=>(1//2, 1//2), "M₂"=>(1//2, 0), "M₃"=>(0, -1//2), "M₄"=>(-1//2, -1//2), "M₅"=>(-1//2, 0), "M₆"=>(0, 1//2)
)
# High-symmetric point in the hexagonal Brillouin zone when the angle between the reciprocal vectors is 60°
const hexagon60°map = Dict(
    "Γ"=>(0, 0), "K"=>(1//3, 1//3), "M"=>(0, 1//2),
    "K₁"=>(1//3, 1//3), "K₂"=>(2//3, -1//3), "K₃"=>(1//3, -2//3), "K₄"=>(-1//3, -1//3), "K₅"=>(-2//3, 1//3), "K₆"=>(-1//3, 2//3),
    "M₁"=>(0, 1//2), "M₂"=>(1//2, 0), "M₃"=>(1//2, -1//2), "M₄"=>(0, -1//2), "M₅"=>(-1//2, 0), "M₆"=>(-1//2, 1//2)
)
"""
    hexagon"P₁-P₂-P₃-..."
    hexagon"P₁-P₂-P₃-..., 120°"
    hexagon"P₁-P₂-P₃-..., 60°"

Construct a tuple of start-stop point pairs for the hexagonal reciprocal space.
"""
macro hexagon_str(str::String)
    str = split(replace(str, " "=>""), ",")
    @assert length(str)∈(1, 2) "@hexagon_str error: wrong pattern."
    length(str)==2 && @assert str[2]∈("120°", "60°") "@hexagon_str error: wrong pattern."
    points = split(str[1], "-")
    @assert length(points)>1 "@hexagon_str error: too few points."
    map = (length(str)==1 || str[2]=="120°") ? hexagon120°map : hexagon60°map
    return (points=ntuple(i->map[points[i]], length(points)), labels=ntuple(i->points[i], length(points)))
end

"""
    ReciprocalCurve{K, S<:SVector} <: ReciprocalSpace{K, S}

A curve in the reciprocal space.
"""
struct ReciprocalCurve{K, S<:SVector} <: ReciprocalSpace{K, S}
    contents::Vector{S}
    function ReciprocalCurve{K}(curve::AbstractVector{<:AbstractVector{<:Number}}) where K
        @assert isa(K, Symbol) "ReciprocalCurve error: K must be a Symbol."
        curve = map(x->SVector{length(x)}(x), curve)
        new{K, eltype(curve)}(curve)
    end
end
@inline VectorSpaceStyle(::Type{<:ReciprocalCurve}) = VectorSpaceEnumerative()

"""
    ReciprocalCurve(curve::AbstractVector{<:NTuple{N, Number}}) where N
    ReciprocalCurve(curve::AbstractVector{<:AbstractVector{<:Number}})
    ReciprocalCurve{K}(curve::AbstractVector{<:NTuple{N, Number}}) where {K, N}
    ReciprocalCurve{K}(curve::AbstractVector{<:AbstractVector{<:Number}}) where K

Construct a curve in the reciprocal space.
"""
@inline ReciprocalCurve(curve::AbstractVector{<:NTuple{N, Number}}) where N = ReciprocalCurve{:k}(curve)
@inline ReciprocalCurve(curve::AbstractVector{<:AbstractVector{<:Number}}) = ReciprocalCurve{:k}(curve)
@inline ReciprocalCurve{K}(curve::AbstractVector{<:NTuple{N, Number}}) where {K, N} = ReciprocalCurve{K}(map(content->SVector{N}(content), curve))

"""
    ReciprocalCurve(reciprocalspace::ReciprocalSpace)

Construct a reciprocal curve from a reciprocal space.
"""
@inline ReciprocalCurve(reciprocalspace::ReciprocalSpace) = ReciprocalCurve{label(reciprocalspace)}(collect(reciprocalspace))

"""
    @recipe plot(curve::ReciprocalCurve)

Define the recipe for the visualization of a reciprocal curve.
"""
@recipe function plot(curve::ReciprocalCurve)
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    legend --> false
    aspect_ratio --> :equal
    minorgrid --> true
    framestyle --> :box
    dim = dimension(curve)
    dim>0 && (xlabel --> label(curve, 1))
    dim>1 && (ylabel --> label(curve, 2))
    dim>2 && (zlabel --> label(curve, 3))
    map(Tuple, curve)
end

# plot utilities
"""
    @recipe plot(reciprocalspace::BrillouinZone, data::AbstractMatrix{<:Number})
    @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractMatrix{<:Number})

Define the recipe for the heatmap visualization of data on a Brillouin/reciprocal zone.
"""
const heatmap = quote
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."
    seriestype --> :heatmap
    titlefontsize --> 10
    aspect_ratio --> :equal
    x, y = range(reciprocalspace, 1), range(reciprocalspace, 2)
    Δx, Δy= x[2]-x[1], y[2]-y[1]
    xlims --> (x[1]-Δx, x[end]+Δx)
    ylims --> (y[1]-Δy, y[end]+Δy)
    clims --> extrema(data)
    xlabel --> label(reciprocalspace, 1)
    ylabel --> label(reciprocalspace, 2)
    x, y, data
end
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractMatrix{<:Number}) = $heatmap
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractMatrix{<:Number}) = $heatmap

"""
    @recipe plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing)

Define the recipe for the scatter visualization of reciprocal points with a series of weights.
"""
@recipe function plot(reciprocalscatter::ReciprocalScatter{K, N}, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing) where {K, N}
    seriestype --> :scatter
    markerstrokewidth --> 0
    legend --> !isnothing(weightlabels)
    point = fractional ? [ntuple(i->0, N)] : [ntuple(i->0, dimension(reciprocalscatter))]
    fractional := fractional
    autolims := false
    for (i, index) in enumerate(axes(weights, 2))
        color = isnothing(weightcolors) ? i : weightcolors[i]
        @series begin
            markersize := atol
            markercolor := color
            label := isnothing(weightlabels) ? "" : weightlabels[i]
            point
        end
        @series begin
            markersize := weights[:, index]*weightmultiplier .+ atol
            markercolor := color
            label := false
            reciprocalscatter
        end
    end
end

"""
    @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number})

Define the recipe for the line visualization of data along a reciprocal path.
"""
const line = quote
    seriestype --> :path
    titlefontsize --> 10
    legend --> false
    xlims --> (0, distance(path))
    minorgrid --> true
    xminorticks --> 10
    yminorticks --> 10
    xticks --> ticks(path)
    xlabel --> label(path)
    [distance(path, i) for i=1:length(path)], data
end
@eval @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}) = $line

"""
    @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing)

Define the recipe for the scatter visualization of data along a reciprocal path with a series of weights.
"""
const scatter = quote
    seriestype --> :scatter
    markerstrokewidth := 0
    legend --> !isnothing(weightlabels)
    xlims --> (0, distance(path))
    alpha --> 0.8
    for i in 1:size(weights, 3)
        @series begin
            markersize := atol
            markercolor := isnothing(weightcolors) ? i : weightcolors[i]
            label := isnothing(weightlabels) ? i : weightlabels[i]
            [(0, 0)]
        end
    end
    for (i, index) in enumerate(axes(weights, 3))
        @series begin
            markersize := weights[:, :, index]*weightmultiplier .+ atol
            markercolor := isnothing(weightcolors) ? i : weightcolors[i]
            label := false
            path, data
        end
    end
end
@eval @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing) = $scatter

"""
    @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number})

Define the recipe for the heatmap visualization of data on the x-y plain with the x axis being a reciprocal path.
"""
@recipe function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number})
    seriestype --> :heatmap
    titlefontsize --> 10
    xticks --> ticks(path)
    xlabel --> label(path)
    xlims --> (0, distance(path))
    ylims --> (minimum(y), maximum(y))
    clims --> extrema(data)
    [distance(path, i) for i=1:length(path)], y, data
end

"""
    @recipe plot(reciprocalspace::BrillouinZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)
    @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)
    @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)

Define the recipe for the heatmap visualization of a series of data on 
1) a Brillouin zone,
2) a reciprocal zone,
3) the x-y plain with the x axis being a reciprocal path.
"""
setup(expr::Expr) = quote
    seriestype := :heatmap
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data)[3])))
    isnothing(ncol) && (ncol = ceil(Int, size(data)[3]/nrow))
    layout := @layout [(nrow, ncol); b{0.05h}]
    colorbar := false
    for i = 1:size(data)[3]
        @series begin
            isnothing(subtitles) || begin
                title := subtitles[i]
                titlefontsize := subtitlefontsize
            end
            subplot := i
            clims := clims
            $expr
        end
    end
    subplot := nrow*ncol+1
    xlims := (minimum(clims), maximum(clims))
    xlabel := ""
    ylims := (0, 1)
    yticks := (0:1, ("", ""))
    ylabel := ""
    LinRange(clims..., 100), [0, 1], [LinRange(clims..., 100)'; LinRange(clims..., 100)']
end
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))
@eval @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(path, y, data[:, :, i])))

# save utilities
"""
    savedlm(filename::AbstractString, x::AbstractVector{<:Number}, y::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, delim='\t')
    savedlm(filename::AbstractString, x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, z::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t')
    savedlm(filename::AbstractString, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t'; fractional::Bool=true)
    savedlm(filename::AbstractString, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}, delim='\t'; fractional::Bool=true)
    savedlm(filename::AbstractString, path::ReciprocalPath, data:AbstractMatrix{<:Number}, delim='\t'; distance::Bool=true)
    savedlm(filename::AbstractString, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}, delim='\t'; distance::Bool=true)
    savedlm(filename::AbstractString, path::ReciprocalPath, y::AbstractVector{<:Number}, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t'; distance::Bool=true)

Save data to delimited files.
"""
function savedlm(filename::AbstractString, x::AbstractVector{<:Number}, y::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, delim='\t')
    @assert length(x)==size(y, 1) "savedlm error: mismatched size of x and y."
    open(filename, "w") do f
        writedlm(f, [x y], delim)
    end
end
function savedlm(filename::AbstractString, x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, z::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t')
    @assert size(z)[1:2]==(length(y), length(x)) "savedlm error: mismatched size of x, y and z."
    open(filename, "w") do f
        new_x = kron(x, ones(length(y)))
        new_y = kron(ones(length(x)), y)
        new_z = reshape(z, length(x)*length(y), :)
        writedlm(f, [new_x new_y new_z], delim)
    end
end
function savedlm(filename::AbstractString, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t'; fractional::Bool=true)
    @assert map(length, shape(reciprocalspace))==reverse(size(data)[1:2]) "savedlm error: mismatched size of reciprocal space and data."
    open(filename, "w") do f
        x = fractional ? matrix(fractionals(reciprocalspace)) : matrix(reciprocalspace)
        y = reshape(data, length(reciprocalspace), :)
        writedlm(f, [x y], delim)
    end
end
function savedlm(filename::AbstractString, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}, delim='\t'; fractional::Bool=true)
    @assert size(weights, 1)==length(reciprocalscatter) "savedlm error: mismatched size of reciprocal scatter points and weights."
    open(filename, "w") do f
        x = fractional ? matrix(fractionals(reciprocalscatter)) : matrix(reciprocalscatter)
        writedlm(f, [x weights], delim)
    end
end
function savedlm(filename::AbstractString, path::ReciprocalPath, data::AbstractMatrix{<:Number}, delim='\t'; distance::Bool=true)
    @assert length(path)==size(data, 1) "savedlm error: mismatched size of path and data."
    open(filename, "w") do f
        x = distance ? [Spatials.distance(path, i) for i=1:length(path), _=1:1] : matrix(path)
        writedlm(f, [x data], delim)
    end
end
function savedlm(filename::AbstractString, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}, delim='\t'; distance::Bool=true)
    @assert length(path)==size(data, 1)==size(weights, 1) && size(data, 2)==size(weights, 2) "savedlm error: mismatched size of path, data and weights."
    open(filename, "w") do f
        x = distance ? kron(ones(size(data, 2)), [Spatials.distance(path, i) for i=1:length(path), _=1:1]) : matrix(kron(ones(size(data, 2)), path))
        y = reshape(data, :)
        z = reshape(weights, length(y), :)
        writedlm(f, [x y z], delim)
    end
end
function savedlm(filename::AbstractString, path::ReciprocalPath, y::AbstractVector{<:Number}, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, delim='\t'; distance::Bool=true)
    @assert (length(y), length(path))==size(data)[1:2] "savedlm error: mismatched size of path, y and data."
    open(filename, "w") do f
        x = distance ? kron([Spatials.distance(path, i) for i=1:length(path), _=1:1], ones(size(data, 1))) : matrix(kron(path, ones(size(data, 1))))
        y = kron(ones(length(path)), y)
        z = reshape(data, prod(size(data)[1:2]), :)
        writedlm(f, [x y z], delim)
    end
end

"""
    matrix(vs::AbstractVector{<:AbstractVector}) -> Matrix

Convert a vector of vector to a matrix.
"""
function matrix(vs::AbstractVector{<:AbstractVector})
    result = zeros(eltype(eltype(vs)), length(vs), length(vs[1]))
    for (i, v) in enumerate(vs)
        result[i, :] = v
    end
    return result
end

end #module
