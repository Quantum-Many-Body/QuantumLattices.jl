module Spatials

using Base.Iterators: flatten, product
using DelimitedFiles: writedlm
using LinearAlgebra: cross, dot, norm
using NearestNeighbors: KDTree, inrange, knn
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series, @layout
using StaticArrays: SVector
using ..QuantumNumbers: Momenta, Momentum, periods
using ..Toolkit: atol, rtol, efficientoperations, CompositeDict, Float, SimpleNamedVectorSpace, Segment, VectorSpaceCartesian, VectorSpaceDirectSummed, VectorSpaceStyle, getcontent

import StaticArrays: SArray
import ..QuantumLattices: decompose, dimension, dtype, expand, kind, reset!
import ..QuantumNumbers: Momentum₁, Momentum₂, Momentum₃
import ..Toolkit: contentnames, shape

export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, bonds, bonds!, icoordinate, isintracell, nneighbor, rcoordinate, save, selectpath, shrink, ticks, xaxis, yaxis, zaxis
export hexagon120°map, hexagon60°map, linemap, rectanglemap, @hexagon_str, @line_str, @rectangle_str

"""
    distance(p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}) -> Number

Get the distance between two points.

!!! note
    Compared to `norm(p₁-p₂)`, this function avoids the memory allocation for `p₁-p₂`, thus is more efficient.
"""
function distance(p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number})
    @assert length(p₁)==length(p₂) "distance error: mismatched length of input vectors."
    result = zero(promote_type(eltype(p₁), eltype(p₂)))
    for i = 1:length(p₁)
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
    result = acos(v[1]/(length(v)==3 ? sqrt(v[1]^2+v[2]^2) : norm(v)))
    (length(v)>1 && v[2]<0) && (result = 2*convert(typeof(result), pi) - result)
    return result
end

"""
    azimuthd(v::AbstractVector{<:Number}) -> Number

Get the azimuth angle in degrees of a vector.
"""
function azimuthd(v::AbstractVector{<:Number})
    @assert length(v)∈(1, 2, 3) "azimuthd error: wrong dimensioned input vector."
    result = acosd(v[1]/(length(v)==3 ? sqrt(v[1]^2+v[2]^2) : norm(v)))
    (length(v)>1 && v[2]<0) && (result = 360 - result)
    return result
end

"""
    polar(v::AbstractVector{<:Number}) -> Number

Get the polar angle in radians of a vector.
"""
function polar(v::AbstractVector{<:Number})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acos(v[3]/norm(v))
end

"""
    polard(v::AbstractVector{<:Number}) -> Number

Get the polar angle in degrees of a vector.
"""
function polard(v::AbstractVector{<:Number})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acosd(v[3]/norm(v))
end

"""
    volume(vectors::AbstractVector{<:SVector}) -> Number
    volume(v::AbstractVector{<:Number}) -> Number
    volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}) -> Number
    volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number}) -> Number

Get the volume spanned by the input vectors.
"""
function volume(vectors::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(vectors)∈(1, 2, 3) "volume error: unsupported number of vectors."
    length(vectors)==1 && return volume(first(vectors))
    length(vectors)==2 && return volume(vectors[1], vectors[2])
    return volume(vectors[1], vectors[2], vectors[3])
end
@inline volume(v::AbstractVector{<:Number}) = norm(v)
function volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number})
    @assert length(v₁)==length(v₂) "volume error: mismatched dimension of vectors."
    @assert length(v₁)∈(1, 2, 3) "volume error: unsupported dimension of vectors."
    length(v₁)==1 && return zero(eltype(v₁))
    length(v₂)==2 && return abs(v₁[1]*v₂[2]-v₁[2]*v₂[1])
    return norm(cross(v₁, v₂))
end
function volume(v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number})
    @assert length(v₁)==length(v₂)==length(v₃) "volume error: mismatched dimension of vectors."
    @assert length(v₁)∈(1, 2, 3) "volume error: unsupported dimension of vectors."
    length(v₁)∈(1, 2) && return zero(eltype(v₁))
    return v₁[1]*(v₂[2]*v₃[3]-v₂[3]*v₃[2]) + v₁[2]*(v₂[3]*v₃[1]-v₂[1]*v₃[3]) + v₁[3]*(v₂[1]*v₃[2]-v₂[2]*v₃[1])
end

"""
    decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}) -> Tuple{Number}
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
    return (sign*n₀/n₁,)
end
function decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number})
    @assert length(v₀)==length(v₁)==length(v₂) "decompose error: mismatched length of input vectors."
    @assert length(v₀)==2 || length(v₀)==3 "decompose error: unsupported dimension($(length(v₀))) of input vectors."
    if length(v₀) == 2
        det = v₁[1]*v₂[2] - v₁[2]*v₂[1]
        x₁ = (v₀[1]*v₂[2]-v₀[2]*v₂[1]) / det
        x₂ = (v₀[2]*v₁[1]-v₀[1]*v₁[2]) / det
    else
        v₃ = SVector{3, promote_type(eltype(v₀), eltype(v₁), eltype(v₂))}(v₁[2]*v₂[3]-v₁[3]*v₂[2], v₁[3]*v₂[1]-v₁[1]*v₂[3], v₁[1]*v₂[2]-v₁[2]*v₂[1])
        x₁, x₂, x₃ = decompose(v₀, v₁, v₂, v₃)
        @assert isapprox(x₃, 0.0, atol=atol, rtol=rtol) "decompose error: insufficient basis vectors."
    end
    return x₁, x₂
end
function decompose(v₀::AbstractVector{<:Number}, v₁::AbstractVector{<:Number}, v₂::AbstractVector{<:Number}, v₃::AbstractVector{<:Number})
    @assert length(v₀)==length(v₁)==length(v₂)==length(v₃) "decompose error: mismatched length of input vectors."
    @assert length(v₀)==3 "decompose error: unsupported dimension($(length(v₀))) of input vectors."
    V = volume(v₁, v₂, v₃)
    r₁ = (v₂[2]*v₃[3]/V-v₂[3]*v₃[2]/V, v₂[3]*v₃[1]/V-v₂[1]*v₃[3]/V, v₂[1]*v₃[2]/V-v₂[2]*v₃[1]/V)
    r₂ = (v₃[2]*v₁[3]/V-v₃[3]*v₁[2]/V, v₃[3]*v₁[1]/V-v₃[1]*v₁[3]/V, v₃[1]*v₁[2]/V-v₃[2]*v₁[1]/V)
    r₃ = (v₁[2]*v₂[3]/V-v₁[3]*v₂[2]/V, v₁[3]*v₂[1]/V-v₁[1]*v₂[3]/V, v₁[1]*v₂[2]/V-v₁[2]*v₂[1]/V)
    return dot(r₁, v₀), dot(r₂, v₀), dot(r₃, v₀)
end
function decompose(v₀::AbstractVector{<:Number}, vs::AbstractVector{<:AbstractVector{<:Number}})
    @assert 0<length(vs)<4 "decompose error: not supported number ($(length(vs))) of vectors."
    result = promote_type(eltype(v₀), eltype(eltype(vs)))[]
    if length(vs)==1
        append!(result, decompose(v₀, vs[1]))
    elseif length(vs)==2
        append!(result, decompose(v₀, vs[1], vs[2]))
    else
        append!(result, decompose(v₀, vs[1], vs[2], vs[3]))
    end
    return result
end

"""
    isintratriangle(
        p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
        vertexes::NTuple{3, Bool}=(true, true, true), edges::NTuple{3, Bool}=(true, true, true), atol::Real=atol, rtol::Real=rtol
    ) -> Bool

Judge whether a point belongs to the interior of a triangle whose vertexes are `p₁`, 'p₂' and `p₃` with the give tolerance. `vertexes` and `edges` define whether the interior should contain the vertexes or edges, respectively.
!!! note
    1. The vertexes are in the order (p₁, p₂, p₃) and the edges are in the order (p1p2, p2p3, p3p1).
    2. The edges do not contain the vertexes.
"""
function isintratriangle(
    p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
    vertexes::NTuple{3, Bool}=(true, true, true), edges::NTuple{3, Bool}=(true, true, true), atol::Real=atol, rtol::Real=rtol
)
    @assert length(p)==length(p₁)==length(p₂)==length(p₃) "isintratriangle error: shape mismatch of input point and triangle."
    @assert length(p)==2 || length(p)==3 "isintratriangle error: unsupported dimension($(length(p))) of input points."
    x = if length(p) == 2
        decompose(SVector(p[1]-p₁[1], p[2]-p₁[2]), SVector(p₂[1]-p₁[1], p₂[2]-p₁[2]), SVector(p₃[1]-p₁[1], p₃[2]-p₁[2]))
    else
        decompose(SVector(p[1]-p₁[1], p[2]-p₁[2], p[3]-p₁[3]), SVector(p₂[1]-p₁[1], p₂[2]-p₁[2], p₂[3]-p₁[3]), SVector(p₃[1]-p₁[1], p₃[2]-p₁[2], p₃[3]-p₁[3]))
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
        result = mapreduce(fapprox, &, decompose(coordinate, vectors[1]))
    elseif length(vectors) == 2
        result = mapreduce(fapprox, &, decompose(coordinate, vectors[1], vectors[2]))
    else
        result = mapreduce(fapprox, &, decompose(coordinate, vectors[1], vectors[2], vectors[3]))
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
    vtype = promote_type(Float, eltype(eltype(vectors)))
    if length(vectors) == 1
        return vector(typeof(vectors), (2*convert(vtype, pi)/mapreduce(vi->vi^2, +, vectors[1])*vectors[1],))
    elseif length(vectors) == 2
        v₁, v₂ = vectors[1], vectors[2]
        @assert length(v₁)==length(v₂)>1 "reciprocals error: mismatched length of input vectors."
        if length(v₁) == 2
            det = 2*convert(vtype, pi) / (v₁[1]*v₂[2]-v₁[2]*v₂[1])
            rv₁ = vector(eltype(vectors), (det*v₂[2], -det*v₂[1]))
            rv₂ = vector(eltype(vectors), (-det*v₁[2], det*v₁[1]))
        else
            v₃ = SVector{3, vtype}(v₁[2]*v₂[3]-v₁[3]*v₂[2], v₁[3]*v₂[1]-v₁[1]*v₂[3], v₁[1]*v₂[2]-v₁[2]*v₂[1])
            V = 2*convert(eltype(v₃), pi) / volume(v₁, v₂, v₃)
            rv₁ = vector(eltype(vectors), (v₂[2]*v₃[3]*V-v₂[3]*v₃[2]*V, v₂[3]*v₃[1]*V-v₂[1]*v₃[3]*V, v₂[1]*v₃[2]*V-v₂[2]*v₃[1]*V))
            rv₂ = vector(eltype(vectors), (v₃[2]*v₁[3]*V-v₃[3]*v₁[2]*V, v₃[3]*v₁[1]*V-v₃[1]*v₁[3]*V, v₃[1]*v₁[2]*V-v₃[2]*v₁[1]*V))
        end
        return vector(typeof(vectors), (rv₁, rv₂))
    elseif length(vectors) == 3
        v₁, v₂, v₃ = vectors
        V = 2*convert(vtype, pi) / volume(v₁, v₂, v₃)
        rv₁ = vector(eltype(vectors), (v₂[2]*v₃[3]*V-v₂[3]*v₃[2]*V, v₂[3]*v₃[1]*V-v₂[1]*v₃[3]*V, v₂[1]*v₃[2]*V-v₂[2]*v₃[1]*V))
        rv₂ = vector(eltype(vectors), (v₃[2]*v₁[3]*V-v₃[3]*v₁[2]*V, v₃[3]*v₁[1]*V-v₃[1]*v₁[3]*V, v₃[1]*v₁[2]*V-v₃[2]*v₁[1]*V))
        rv₃ = vector(eltype(vectors), (v₁[2]*v₂[3]*V-v₁[3]*v₂[2]*V, v₁[3]*v₂[1]*V-v₁[1]*v₂[3]*V, v₁[1]*v₂[2]*V-v₁[2]*v₂[1]*V))
        return vector(typeof(vectors), (rv₁, rv₂, rv₃))
    end
end
@inline vector(::Type{<:Vector}, v::Tuple) = collect(v)
@inline vector(::Type{<:SVector}, v::Tuple) = SVector(v)

"""
    rotate(vector::AbstractVector{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))) -> Vector{<:Number}
    rotate(cluster::AbstractMatrix{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))) -> Matrix{<:Number}

Get a rotated vector/cluster of the original one by a certain angle around an axis.

The axis is determined by a point it gets through (`nothing` can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.
!!! note
    1. The result is given by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
    2. Only 2 and 3 dimensional vectors can be rotated.
    3. When the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.
"""
@inline function rotate(vector::AbstractVector{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0)))
    return reshape(rotate(reshape(vector, :, 1), angle; axis=axis), :)
end
function rotate(cluster::AbstractMatrix{<:Number}, angle::Number; axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0)))
    @assert size(cluster, 1)∈(2, 3) "rotate error: only 2 and 3 dimensional vectors can be rotated."
    datatype = promote_type(eltype(cluster), typeof(angle), Float)
    center, theta, phi = (isnothing(axis[1]) ? zeros(datatype, size(cluster, 1)) : axis[1]), axis[2][1], axis[2][2]
    @assert length(center)==size(cluster, 1) "rotate error: mismatched shape of the input cluster and the point on axis."
    if length(center) == 2
        @assert isapprox(theta, 0, atol=atol, rtol=rtol) && isapprox(phi, 0, atol=atol, rtol=rtol) "rotate error: both the polar and azimuth of the axis for 2d vectors must be 0."
    end
    cosθ, sinθ = cos(angle), sin(angle)
    k, w = SVector{3, datatype}(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)), zeros(datatype, 3)
    result = zeros(datatype, size(cluster))
    for i = 1:size(cluster, 2)
        for j = 1:size(cluster, 1)
            w[j] = cluster[j, i] - center[j]
        end
        inner = dot(k, w)
        outer = (k[2]*w[3]-k[3]*w[2], k[3]*w[1]-k[1]*w[3], k[1]*w[2]-k[2]*w[1])
        for j = 1:size(cluster, 1)
            result[j, i] = w[j]*cosθ + outer[j]*sinθ + k[j]*inner*(1-cosθ) + center[j]
        end
    end
    return result
end

"""
    translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number}) -> Matrix{vector|>eltype}

Get the translated cluster of the original one by a vector.
"""
@inline translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number}) = cluster .+ reshape(vector, (vector|>length, 1))

"""
    tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations) -> Matrix{<:Number}

Tile a supercluster by translations of the input cluster.

Basically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by `vectors` and each set of the translation indices contained in `translations`. When translation vectors are empty, a copy of the original cluster will be returned.
"""
function tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations)
    length(vectors)==0 && return copy(cluster)
    length(translations)>0 && @assert length(vectors)==length(first(translations)) "tile error: mismatched shape of input vectors and translations."
    datatype = promote_type(eltype(cluster), eltype(eltype(vectors)), Float)
    supercluster = zeros(datatype, size(cluster, 1), size(cluster, 2)*length(translations))
    disp = zeros(datatype, size(cluster, 1))
    for (i, translation) in enumerate(translations)
        for i = 1:length(disp)
            disp[i] = zero(datatype)
            for j = 1:length(vectors)
                disp[i] += vectors[j][i] * translation[j]
            end
        end
        for j = 1:size(cluster, 2)
            col = (i-1)*size(cluster, 2) + j
            for row = 1:size(cluster, 1)
                supercluster[row, col] = cluster[row, j] + disp[row]
            end
        end
    end
    return supercluster
end

"""
    minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Int=1; coordination::Int=12) -> Vector{Float}

Use kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.

When the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds across the periodic boundaries, the cluster will be pre-translated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The `coordination` parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain `Inf`. This is a sign that you may need a larger `coordination`. Another situation that `Inf` appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.
"""
function minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Int=1; coordination::Int=12)
    @assert nneighbor>=0 "minimumlengths error: input nneighbor must be non negative."
    result = [Inf for i = 1:nneighbor]
    if size(cluster, 2) > 0
        cluster, vectors = convert(Matrix{Float}, cluster), convert(Vector{Vector{Float}}, vectors)
        translations = reshape(product((-nneighbor:nneighbor for i = 1:length(vectors))...)|>collect, :)
        for translation in translations
            if length(translation)>0 && mapreduce(≠(0), |, translation)
                remove = ntuple(i->-translation[i], length(translation))
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
@inline Neighbors(lengths::Vector{<:Number}, ordinals=0:length(lengths)-1) = Neighbors(ordinal=>length for (ordinal, length) in zip(ordinals, lengths))
@inline Base.max(neighbors::Neighbors) = max(values(neighbors)...)
@inline nneighbor(neighbors::Neighbors{<:Integer}) = max(keys(neighbors)...)
@inline nneighbor(neighbors::Neighbors) = length(neighbors)

"""
    interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Neighbors) -> Vector{Tuple{Int, Int, Int}}

Use kdtree to get the intercluster nearest neighbors.
"""
function interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Neighbors)
    @assert size(cluster₁, 1)==size(cluster₂, 1) "interlinks error: mismatched space dimension of input clusters."
    result = Tuple{Int, Int, Int}[]
    length(neighbors)==0 && return result
    for (i, indices) in enumerate(inrange(KDTree(convert(Matrix{Float}, cluster₂)), convert(Matrix{Float}, cluster₁), max(neighbors)+atol, true))
        for j in indices
            dist = zero(promote_type(eltype(cluster₁), eltype(cluster₂)))
            for k = 1:size(cluster₁, 1)
                dist = dist + (cluster₂[k, j]-cluster₁[k, i])^2
            end
            dist = sqrt(dist)
            for (nb, len) in neighbors
                if isapprox(len, dist, atol=atol, rtol=rtol)
                    push!(result, (nb, j, i))
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
    Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, N)) where N
    Point(site::Integer, rcoordinate::AbstractVector{<:Number}, icoordinate::AbstractVector{<:Number}=zero(SVector{length(rcoordinate), Int}))

Construct a labeled point.
"""
@inline function Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, N)) where N
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
    dtype(point::Point)
    dtype(::Type{<:Point{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a point.
"""
@inline dtype(point::Point) = dtype(typeof(point))
@inline dtype(::Type{<:Point{N, D} where N}) where {D<:Number} = D

"""
    isintracell(point::Point) -> Bool

Judge whether a point is intra the unitcell.
"""
@inline isintracell(point::Point) = isapprox(norm(point.icoordinate), 0.0, atol=atol, rtol=rtol)

"""
    Bond{K, P<:Point}

A generic bond, which could contains several points.
"""
struct Bond{K, P<:Point}
    kind::K
    points::Vector{P}
end
@inline Base.:(==)(bond₁::Bond, bond₂::Bond) = ==(efficientoperations, bond₁, bond₂)
@inline Base.isequal(bond₁::Bond, bond₂::Bond) = isequal(efficientoperations, bond₁, bond₂)
@inline Base.show(io::IO, bond::Bond) = @printf io "Bond(%s, %s)" repr(bond.kind) join(map(string, bond.points), ", ")

"""
    Bond(point::Point)
    Bond(kind, point₁::Point, point₂::Point, points::Point...)

Construct a bond.
"""
@inline Bond(point::Point) = Bond(0, [point])
@inline Bond(kind, point₁::Point, point₂::Point, points::Point...) = Bond(kind, [point₁, point₂, points...])

"""
    dimension(bond::Bond) -> Int
    dimension(::Type{<:Bond{K, P} where K}) where {P<:Point} -> Int

Get the space dimension of a concrete bond.
"""
@inline dimension(bond::Bond) = dimension(typeof(bond))
@inline dimension(::Type{<:Bond{K, P} where K}) where {P<:Point} = dimension(P)

"""
    dtype(bond::Bond)
    dtype(::Type{<:Bond{K, P} where K}) where {P<:Point}

Get the data type of the coordinates of the points contained in a generic bond.
"""
@inline dtype(bond::Bond) = dtype(typeof(bond))
@inline dtype(::Type{<:Bond{K, P} where K}) where {P<:Point} = dtype(P)

"""
    length(bond::Bond) -> Int

Get the number of points contained in a generic bond.
"""
@inline Base.length(bond::Bond) = length(bond.points)

"""
    eltype(bond::Bond)
    eltype(::Type{<:Bond{K, P} where K}) where {P<:Point}

Get the point type contained in a generic bond.
"""
@inline Base.eltype(bond::Bond) = eltype(typeof(bond))
@inline Base.eltype(::Type{<:Bond{K, P} where K}) where {P<:Point} = P

"""
    iterate(bond::Bond, state=1)

Iterate over the points contained in a generic bond.
"""
@inline Base.iterate(bond::Bond, state=1) = state<=length(bond) ? (bond[state], state+1) : nothing

"""
    reverse(bond::Bond) -> Bond

Get the reversed bond.
"""
@inline Base.reverse(bond::Bond) = Bond(bond.kind, reverse(bond.points))

"""
    getindex(bond::Bond, i::Integer) -> Point

Get the ith point contained in a generic bond.
"""
@inline Base.getindex(bond::Bond, i::Integer) = bond.points[i]

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
    dtype(lattice::AbstractLattice)
    dtype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a lattice.
"""
@inline dtype(lattice::AbstractLattice) = dtype(typeof(lattice))
@inline dtype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number} = D

"""
    length(lattice::AbstractLattice) -> Int

Get the number of points contained in a lattice.
"""
@inline Base.length(lattice::AbstractLattice) = size(getcontent(lattice, :coordinates))[2]

"""
    getindex(lattice::AbstractLattice, i::Integer) -> SVector

Get the ith coordinate.
"""
@inline Base.getindex(lattice::AbstractLattice, i::Integer) = SVector{dimension(lattice), dtype(lattice)}(ntuple(j->lattice.coordinates[j, i], Val(dimension(lattice))))

"""
    reciprocals(lattice::AbstractLattice) -> Vector{<:SVector}

Get the reciprocal translation vectors of the dual lattice.
"""
@inline reciprocals(lattice::AbstractLattice) = reciprocals(getcontent(lattice, :vectors))

"""
    Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Int=12)

Get the neighbor vs. bond length map of a lattice up to the `nneighbor`th order.
"""
@inline Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Int=12) = Neighbors(minimumlengths(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), nneighbor; coordination=coordination))

"""
    bonds(lattice::AbstractLattice, nneighbor::Int; coordination::Int=12) -> Vector{Bond{Int, Point{dimension(lattice), dtype(lattice)}}}
    bonds(lattice::AbstractLattice, neighbors::Neighbors) -> Vector{Bond{keytype(neighbors), Point{dimension(lattice), dtype(lattice)}}}

Get the required bonds of a lattice.
"""
@inline bonds(lattice::AbstractLattice, nneighbor::Int; coordination::Int=12) = bonds!(Bond{Int, Point{dimension(lattice), dtype(lattice)}}[], lattice, nneighbor; coordination=coordination)
@inline bonds(lattice::AbstractLattice, neighbors::Neighbors) = bonds!(Bond{keytype(neighbors), Point{dimension(lattice), dtype(lattice)}}[], lattice, neighbors)

"""
    bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Int; coordination::Int=12)
    bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors) -> typeof(bonds)

Get the required bonds of a lattice and append them to the input bonds.
"""
@inline bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Int; coordination::Int=12) = bonds!(bonds, lattice, Neighbors(lattice, nneighbor; coordination=coordination))
function bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors)
    origin = zero(SVector{dimension(lattice), dtype(lattice)})
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
    translations = reshape(product((-nnb:nnb for i = 1:length(getcontent(lattice, :vectors)))...)|>collect, :)
    for translation in translations
        remove = ntuple(i->-translation[i], length(translation))
        filter!(≠(remove), translations)
    end
    if length(translations) > 0
        superrcoordinates = tile(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), translations)
        supericoordinates = tile(zero(getcontent(lattice, :coordinates)), getcontent(lattice, :vectors), translations)
        for (k, index₁, index₂) in interlinks(getcontent(lattice, :coordinates), superrcoordinates, neighbors)
            rcoordinate = SVector{dimension(lattice), dtype(lattice)}(ntuple(j->superrcoordinates[j, index₁], Val(dimension(lattice))))
            icoordinate = SVector{dimension(lattice), dtype(lattice)}(ntuple(j->supericoordinates[j, index₁], Val(dimension(lattice))))
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
    (isa(neighbors, Int) || 0∈keys(neighbors)) && begin
        @series begin
            seriestype := :scatter
            coordinates = NTuple{dimension(lattice), dtype(lattice)}[]
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
    @series begin
        data = Vector{Float64}[]
        for i = 1:dimension(lattice)
            push!(data, Float64[])
        end
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
function Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    vectors = isnothing(vectors) ? SVector{0, SVector{N, eltype(eltype(coordinates))}}() : vectorconvert(vectors)
    coordinates = [coordinates[j][i] for i=1:N, j=1:length(coordinates)]
    return Lattice(name, coordinates, vectors)
end
function Lattice(coordinates::AbstractVector{<:Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing)
    coordinates = hcat(coordinates...)
    vectors = isnothing(vectors) ? SVector{0, SVector{size(coordinates)[1], eltype(coordinates)}}() : vectorconvert(vectors)
    return Lattice(name, coordinates, vectors)
end
@inline vectorconvert(vectors::SVector{N, <:SVector}) where N = vectors
@inline vectorconvert(vectors::AbstractVector{<:SVector}) = convert(SVector{length(vectors), SVector{length(eltype(vectors)), eltype(eltype(vectors))}}, vectors)
@inline vectorconvert(vectors::AbstractVector{<:AbstractVector}) = convert(SVector{length(vectors), SVector{length(first(vectors)), eltype(eltype(vectors))}}, vectors)

"""
    Lattice(lattice::AbstractLattice, ranges::NTuple{N, Int}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N)); mode::Symbol=:nonnegative) where N
    Lattice(lattice::AbstractLattice, ranges::NTuple{N, UnitRange{Int}}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N))) where N

Construct a lattice from the translations of another.
"""
function Lattice(lattice::AbstractLattice, ranges::NTuple{N, Int}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N)); mode::Symbol=:nonnegative) where N
    @assert mode∈(:center, :nonnegative) "Lattice error: wrong mode($(repr(mode)))."
    return Lattice(lattice, mode==:center ? map(i->-floor(Int, (i-1)/2):-floor(Int, (i-1)/2)+i-1, ranges) : map(i->0:i-1, ranges), boundaries)
end
function Lattice(lattice::AbstractLattice, ranges::NTuple{N, UnitRange{Int}}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N))) where N
    @assert all(map(in(('P', 'O', 'p', 'o')), boundaries)) "Lattice error: boundary conditions must be either 'P'/'p' for 'periodic' or 'O'/'o' for 'open'."
    @assert length(boundaries)==N "Lattice error: mismatched number of ranges and boundaries conditions."
    boundaries = map(uppercase, boundaries)
    name = Symbol(@sprintf "%s%s" getcontent(lattice, :name) join([@sprintf("%s%s%s", boundary=='P' ? "[" : "(", range, boundary=='P' ? "]" : ")") for (range, boundary) in zip(ranges, boundaries)]))
    coordinates = tile(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), product(ranges...))
    vectors = SVector{dimension(lattice), dtype(lattice)}[]
    for (i, vector) in enumerate(getcontent(lattice, :vectors))
        boundaries[i]=='P' && push!(vectors, vector*length(ranges[i]))
    end
    return Lattice(name, coordinates, vectorconvert(vectors))
end

"""
    expand(momentum::Momentum, reciprocals::AbstractVector{<:AbstractVector}) -> eltype(reciprocals)

Expand the momentum from integral values to real values with the given reciprocals.
"""
@inline function expand(momentum::Momentum, reciprocals::AbstractVector{<:AbstractVector})
    p = periods(momentum)
    @assert length(p)==length(reciprocals) "expand error: mismatched momentum and reciprocals."
    result = zero(first(reciprocals))
    for i = 1:length(p)
        result += reciprocals[i]*(momentum[i]//p[i])
    end
    return result
end

"""
    Momentum₁{N}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where N
    Momentum₂{N₁, N₂}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where {N₁, N₂}
    Momentum₃{N₁, N₂, N₃}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where {N₁, N₂, N₃}

Construct a quantum momentum by the coordinates.
"""
function Momentum₁{N}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where N
    @assert length(reciprocals)==1 "Momentum₁ error: mismatched length of reciprocals."
    k = decompose(momentum, reciprocals[1])[1]*N
    i = round(Int, k)
    @assert isapprox(i, k; atol=atol, rtol=rtol) "Momentum₁ error: input momentum not on grid."
    return Momentum₁{N}(i)
end
function Momentum₂{N₁, N₂}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where {N₁, N₂}
    @assert length(reciprocals)==2 "Momentum₂ error: mismatched length of reciprocals."
    k₁, k₂ = decompose(momentum, reciprocals[1], reciprocals[2])
    k₁, k₂ = k₁*N₁, k₂*N₂
    i₁, i₂ = round(Int, k₁), round(Int, k₂)
    @assert isapprox(i₁, k₁; atol=atol, rtol=rtol) && isapprox(i₂, k₂; atol=atol, rtol=rtol) "Momentum₂ error: input momentum not on grid."
    return Momentum₂{N₁, N₂}(i₁, i₂)
end
function Momentum₃{N₁, N₂, N₃}(momentum::AbstractVector, reciprocals::AbstractVector{<:AbstractVector}) where {N₁, N₂, N₃}
    @assert length(reciprocals)==3 "Momentum₃ error: mismatched length of reciprocals."
    k₁, k₂, k₃ = decompose(momentum, reciprocals[1], reciprocals[2], reciprocals[3])
    k₁, k₂, k₃ = k₁*N₁, k₂*N₂, k₃*N₃
    i₁, i₂, i₃ = round(Int, k₁), round(Int, k₂), round(Int, k₃)
    @assert isapprox(i₁, k₁; atol=atol, rtol=rtol) && isapprox(i₂, k₂; atol=atol, rtol=rtol) && isapprox(i₃, k₃; atol=atol, rtol=rtol) "Momentum₃ error: input momentum not on grid."
    return Momentum₃{N₁, N₂, N₃}(i₁, i₂, i₃)
end

"""
    ReciprocalSpace{K, P<:SVector} <: SimpleNamedVectorSpace{K, P}

Abstract type of reciprocal spaces.
"""
abstract type ReciprocalSpace{K, P<:SVector} <: SimpleNamedVectorSpace{K, P} end
@inline dtype(reciprocalspace::ReciprocalSpace) = dtype(typeof(reciprocalspace))
@inline dtype(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector} = eltype(P)
@inline dimension(reciprocalspace::ReciprocalSpace) = dimension(typeof(reciprocalspace))
@inline dimension(::Type{<:ReciprocalSpace{K, P} where K}) where {P<:SVector} = length(P)

"""
    @recipe plot(reciprocalspace::ReciprocalSpace)

Define the recipe for the visualization of a reciprocal space.
"""
@recipe function plot(reciprocalspace::ReciprocalSpace)
    title --> string(nameof(typeof(reciprocalspace)))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    seriestype := :scatter
    map(Tuple, reciprocalspace)
end

"""
    BrillouinZone{K, P<:Momentum, S<:SVector, N} <: ReciprocalSpace{K, S}

The Brillouin zone of a lattice.
"""
struct BrillouinZone{K, P<:Momentum, S<:SVector, N} <: ReciprocalSpace{K, S}
    reciprocals::SVector{N, S}
    function BrillouinZone{K, P}(reciprocals::SVector{N, <:SVector}) where {K, P<:Momentum, N}
        @assert isa(K, Symbol) "BrillouinZone error: K must be a Symbol."
        @assert fieldcount(P)==N "BrillouinZone error: mismatched momentum and reciprocals."
        new{K, P, eltype(reciprocals), N}(reciprocals)
    end
end
@inline Base.:(==)(bz₁::BrillouinZone, bz₂::BrillouinZone) = keytype(bz₁)==keytype(bz₂) && bz₁.reciprocals==bz₂.reciprocals
@inline Base.isequal(bz₁::BrillouinZone, bz₂::BrillouinZone) = isequal(keytype(bz₁), keytype(bz₂)) && isequal(bz₁.reciprocals, bz₂.reciprocals)
@inline Base.hash(brillouinzone::BrillouinZone, h::UInt) = hash((ntuple(i->brillouinzone.reciprocals[i], brillouinzone|>keytype|>fieldcount|>Val), brillouinzone|>keytype|>periods), h)
@inline VectorSpaceStyle(::Type{<:BrillouinZone}) = VectorSpaceCartesian()
@inline shape(::BrillouinZone{K, P}) where {K, P<:Momentum} = shape(Momenta(P))
@inline SArray(index::CartesianIndex, brillouinzone::BrillouinZone{K, P}) where {K, P<:Momentum} = expand(P(index, Momenta(P)), brillouinzone.reciprocals)
@inline Base.keys(::BrillouinZone{K, P}) where {K, P<:Momentum} = Momenta(P)
@inline Base.keytype(brillouinzone::BrillouinZone) = keytype(typeof(brillouinzone))
@inline Base.keytype(::Type{<:BrillouinZone{K, P} where K}) where {P<:Momentum} = P
@inline xaxis(brillouinzone::BrillouinZone) = (n=periods(keytype(brillouinzone))[1]; collect(Float64, 0:(n-1))/n)
@inline yaxis(brillouinzone::BrillouinZone) = (n=periods(keytype(brillouinzone))[2]; collect(Float64, 0:(n-1))/n)
@inline zaxis(brillouinzone::BrillouinZone) = (n=periods(keytype(brillouinzone))[3]; collect(Float64, 0:(n-1))/n)

"""
    BrillouinZone(reciprocals::AbstractVector{<:AbstractVector}, nk)
    BrillouinZone{K}(reciprocals::AbstractVector{<:AbstractVector}, nk) where K
    BrillouinZone(::Type{P}, reciprocals::AbstractVector{<:AbstractVector}) where {P<:Momentum}
    BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector{<:AbstractVector}) where {K, P<:Momentum}

Construct a Brillouin zone.
"""
@inline BrillouinZone(reciprocals::AbstractVector{<:AbstractVector}, nk) = BrillouinZone{:k}(reciprocals, nk)
@inline BrillouinZone(::Type{P}, reciprocals::AbstractVector{<:AbstractVector}) where {P<:Momentum} = BrillouinZone{:k}(P, reciprocals)
@inline function BrillouinZone{K}(reciprocals::AbstractVector{<:AbstractVector}, nk) where K
    isa(nk, Integer) && (nk = ntuple(i->nk, length(reciprocals)))
    @assert length(nk)==length(reciprocals) "BrillouinZone error: mismatched number of reciprocals and periods of momenta."
    length(reciprocals)==1 && return BrillouinZone{K}(Momentum₁{nk[1]}, reciprocals)
    length(reciprocals)==2 && return BrillouinZone{K}(Momentum₂{nk[1], nk[2]}, reciprocals)
    length(reciprocals)==3 && return BrillouinZone{K}(Momentum₃{nk[1], nk[2], nk[3]}, reciprocals)
    error("BrillouinZone error: only 1d, 2d and 3d are supported.")
end
@inline BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector{<:AbstractVector}) where {K, P<:Momentum} = BrillouinZone{K, P}(vectorconvert(reciprocals))

"""
    ReciprocalZone{K, S<:SVector, V<:Number} <: ReciprocalSpace{K, S}

A zone in the reciprocal space.
"""
struct ReciprocalZone{K, S<:SVector, N, V<:Number} <: ReciprocalSpace{K, S}
    reciprocals::SVector{N, S}
    bounds::NTuple{N, Segment{V}}
    volume::V
    function ReciprocalZone{K}(reciprocals::SVector{N, <:SVector}, bounds::NTuple{N, <:Segment}) where {K, N}
        @assert isa(K, Symbol) "ReciprocalZone error: K must be a Symbol."
        ratio = one(promote_type(eltype(eltype(reciprocals)), eltype(eltype(bounds))))
        for i = 1:length(bounds)
            ratio = ratio*(bounds[i].stop-bounds[i].start)
        end
        v = ratio*volume(reciprocals)
        new{K, eltype(reciprocals), N, typeof(v)}(reciprocals, bounds, v)
    end
end
@inline VectorSpaceStyle(::Type{<:ReciprocalZone}) = VectorSpaceCartesian()
@inline shape(reciprocalzone::ReciprocalZone) = map(bound->1:length(bound), reverse(reciprocalzone.bounds))
@inline function SArray(index::CartesianIndex, reciprocalzone::ReciprocalZone)
    result = zero(eltype(reciprocalzone))
    for (reciprocal, bound, i) in zip(reciprocalzone.reciprocals, reciprocalzone.bounds, reverse(index.I))
        result += reciprocal*bound[i]
    end
    return result
end
@inline xaxis(reciprocalzone::ReciprocalZone) = collect(Float64, reciprocalzone.bounds[1])
@inline yaxis(reciprocalzone::ReciprocalZone) = collect(Float64, reciprocalzone.bounds[2])
@inline zaxis(reciprocalzone::ReciprocalZone) = collect(Float64, reciprocalzone.bounds[3])

"""
    shrink(reciprocalzone::ReciprocalZone{K}, ranges::Vararg{OrdinalRange{<:Integer}, N}) where {K, N} -> ReciprocalZone

Shrink a reciprocal zone.
"""
function shrink(reciprocalzone::ReciprocalZone{K}, ranges::Vararg{OrdinalRange{<:Integer}, N}) where {K, N}
    @assert length(ranges)==length(reciprocalzone.reciprocals) "shrink error: mismatched number of ranges and reciprocals."
    return ReciprocalZone{K}(reciprocalzone.reciprocals, ntuple(i->reciprocalzone.bounds[i][ranges[i]], N|>Val))
end

"""
    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}; length=100, ends=(true, false))
    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false))
    ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}, bounds::Union{Tuple{Vararg{Pair{<:Number, <:Number}}}, Vector{<:Pair{<:Number, <:Number}}}; length=100, ends=(true, false))

    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}; length=100, ends=(true, false)) where K
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false)) where K
    ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}, bounds::Union{Tuple{Vararg{Pair{<:Number, <:Number}}}, Vector{<:Pair{<:Number, <:Number}}}; length=100, ends=(true, false)) where K

Construct a rectangular zone in the reciprocal space.
"""
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals); length=length, ends=ends)
end
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}, bound::Pair{<:Number, <:Number}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals), bound, bounds...; length=length, ends=ends)
end
@inline function ReciprocalZone(reciprocals::AbstractVector{<:AbstractVector}, bounds::Union{Tuple{Vararg{Pair{<:Number, <:Number}}}, Vector{<:Pair{<:Number, <:Number}}}; length=100, ends=(true, false))
    return ReciprocalZone{:k}(vectorconvert(reciprocals), bounds; length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}; length=100, ends=(true, false)) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), ntuple(i->-1//2=>1//2, Base.length(reciprocals)); length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}, bound::Pair{<:Number, <:Number}, bounds::Pair{<:Number, <:Number}...; length=100, ends=(true, false)) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), (bound, bounds...); length=length, ends=ends)
end
@inline function ReciprocalZone{K}(reciprocals::AbstractVector{<:AbstractVector}, bounds::Union{Tuple{Vararg{Pair{<:Number, <:Number}}}, Vector{<:Pair{<:Number, <:Number}}}; length=100, ends=(true, false)) where K
    @assert Base.length(reciprocals)==Base.length(bounds) "ReciprocalZone error: the numbers of the input reciprocals($(Base.length(reciprocals))) and bounds($(Base.length(bounds))) do not match."
    isa(length, Integer) && (length = map(i->length, 1:Base.length(bounds)))
    isa(ends, NTuple{2, Bool}) && (ends = map(i->ends, 1:Base.length(bounds)))
    @assert Base.length(bounds)==Base.length(length) "ReciprocalZone error: the number of length should be ($(Base.length(bounds))) if it is not `Integer`."
    @assert Base.length(bounds)==Base.length(ends) "ReciprocalZone error: the number of ends should be ($(Base.length(bounds))) if it is not `NTuple{2, Bool}`."
    segments = ntuple(i->Segment(bounds[i].first, bounds[i].second, length[i]; ends=ends[i]), Base.length(bounds))
    return ReciprocalZone{K}(vectorconvert(reciprocals), segments)
end

"""
    ReciprocalZone(brillouinzone::BrillouinZone)

Construct a reciprocal zone from a Brillouin zone.
"""
@inline function ReciprocalZone(brillouinzone::BrillouinZone{K, P}) where {K, P<:Momentum}
    return ReciprocalZone{K}(brillouinzone.reciprocals, map(length->Segment(0, 1, length; ends=(true, false)), periods(P)))
end

"""
    ReciprocalPath{K, S<:SVector, N, R} <: ReciprocalSpace{K, S}

A path in the reciprocal space.
"""
struct ReciprocalPath{K, S<:SVector, N, R} <: ReciprocalSpace{K, S}
    contents::NTuple{N, Vector{S}}
    labels::NTuple{N, Pair{R, R}}
    function ReciprocalPath{K}(contents::NTuple{N, AbstractVector{<:AbstractVector}}, labels::NTuple{N, Pair{R, R}}) where {K, N, R}
        @assert isa(K, Symbol) "ReciprocalPath error: K must be a Symbol."
        new{K, eltype(eltype(contents)), N, R}(map(vectorconvert, contents), labels)
    end
end
@inline ReciprocalPath(contents::Tuple{Vararg{AbstractVector{<:AbstractVector}}}, labels::Tuple{Vararg{Pair}}) = ReciprocalPath{:k}(contents, labels)
@inline contentnames(::Type{<:ReciprocalPath}) = (:contents, :labels)
@inline VectorSpaceStyle(::Type{<:ReciprocalPath}) = VectorSpaceDirectSummed()

points2segments(points::NTuple{M}) where M = ntuple(i->points[i]=>points[i+1], Val(M-1))
"""
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing)
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, points::NTuple{N, Number}...; labels=points, length=100, ends=nothing) where N
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, length=100, ends=nothing) where N
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, length=100, ends=nothing) where N
    ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, length=100, ends=nothing) where N

    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing) where K
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, points::NTuple{N, Number}...; labels=points, length=100, ends=nothing) where {K, N}
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, length=100, ends=nothing) where {K, N}
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, length=100, ends=nothing) where {K, N}
    ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, length=100, ends=nothing) where {K, N}

Construct a path in the reciprocal space.

When `length` is an integer, it specifies the length of each segment except for the last whose length will be `length+1`.
When `ends` is `nothing`, the start point will be included while the end point will be not for each segment except for the last whose both points will be included.
"""
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing)
    return ReciprocalPath(reciprocals, contents.points...; labels=contents.labels, length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, points::NTuple{N, Number}...; labels=points, length=100, ends=nothing) where N
    return ReciprocalPath(reciprocals, points; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, length=100, ends=nothing) where N
    return ReciprocalPath(reciprocals, points2segments(points)...; labels=points2segments(Tuple(labels)), length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, length=100, ends=nothing) where N
    return ReciprocalPath(reciprocals, segments; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath(reciprocals::AbstractVector{<:AbstractVector}, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, length=100, ends=nothing) where N
    return ReciprocalPath{:k}(reciprocals, segments...; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; length=100, ends=nothing) where K
    return ReciprocalPath{K}(reciprocals, contents.points...; labels=contents.labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, points::NTuple{N, Number}...; labels=points, length=100, ends=nothing) where {K, N}
    return ReciprocalPath{K}(reciprocals, points; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, length=100, ends=nothing) where {K, N}
    return ReciprocalPath{K}(reciprocals, points2segments(points)...; labels=points2segments(Tuple(labels)), length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, length=100, ends=nothing) where {K, N}
    return ReciprocalPath{K}(reciprocals, segments; labels=labels, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector{<:AbstractVector}, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, length=100, ends=nothing) where {K, N}
    @assert Base.length(reciprocals)==N "ReciprocalPath error: mismatched input reciprocals and segments."
    isa(length, Integer) && (length = ntuple(i->i<fieldcount(typeof(segments)) ? (segments[i].second==segments[i+1].first ? length : length+1) : length+1, fieldcount(typeof(segments))))
    isnothing(ends) && (ends = ntuple(i->i<fieldcount(typeof(segments)) ? (true, segments[i].second==segments[i+1].first ? false : true) : (true, true), fieldcount(typeof(segments))))
    @assert fieldcount(typeof(segments))==fieldcount(typeof(length)) "ReciprocalPath error: the number of length should be ($(fieldcount(typeof(segments)))) if it is not `Integer`."
    @assert fieldcount(typeof(segments))==fieldcount(typeof(ends)) "ReciprocalPath error: the number of ends should be ($(fieldcount(typeof(segments)))) if it is not `NTuple{2, Bool}`."
    reciprocals = vectorconvert(reciprocals)
    segments = ntuple(i->Segment(segments[i].first, segments[i].second, length[i]; ends=ends[i]), fieldcount(typeof(segments)))
    contents = ntuple(i->zeros(SVector{Base.length(eltype(reciprocals)), promote_type(eltype(eltype(reciprocals)), eltype(eltype(eltype(segments))))}, Base.length(segments[i])), fieldcount(typeof(segments)))
    for (m, segment) in enumerate(segments)
        for (n, index) in enumerate(segment)
            for i = 1:Base.length(index)
                contents[m][n] += reciprocals[i]*index[i]
            end
        end
    end
    return ReciprocalPath{K}(contents, Tuple(labels))
end

"""
    step(path::ReciprocalPath, i::Int) -> dtype(path)

Get the step between the ith and the (i+1)th points in the path.
"""
@inline Base.step(path::ReciprocalPath, i::Int) = distance(path[i], path[i+1])

"""
    distance(path::ReciprocalPath) -> dtype(path)
    distance(path::ReciprocalPath, i::Int) -> dtype(path)
    distance(path::ReciprocalPath, i::Int, j::Int) -> dtype(path)

Get the distance
1) of the total path,
2) from the start point to the ith point in the path,
3) from the ith point to the jth point in the path (when i is greater than j, the value is negative).
"""
@inline distance(path::ReciprocalPath) = distance(path, 1, length(path))
@inline distance(path::ReciprocalPath, i::Int) = distance(path, 1, i)
function distance(path::ReciprocalPath, i::Int, j::Int)
    i==j && return zero(dtype(path))
    i>j && return -distance(path, j, i)
    dimsum = cumsum(map(length, path.contents))
    i₁ = searchsortedfirst(dimsum, i)
    j₁ = searchsortedfirst(dimsum, j)
    i₂ = i₁>1 ? (i-dimsum[i₁-1]) : i
    j₂ = j₁>1 ? (j-dimsum[j₁-1]) : j
    i₁==j₁ && return distance(path.contents[i₁][j₁], path.contents[i₂][j₂])
    result = distance(path.contents[i₁][i₂], path.contents[i₁][end]) + distance(path.contents[j₁][1], path.contents[j₁][j₂])
    for l = (i₁+1):(j₁-1)
        path.labels[l-1].second==path.labels[l].first && (result += distance(path.contents[l-1][end], path.contents[l][1]))
        result += distance(path.contents[l][1], path.contents[l][end])
    end
    path.labels[j₁-1].second==path.labels[j₁].first && (result += distance(path.contents[j₁-1][end], path.contents[j₁][1]))
    return result
end

"""
    ticks(path::ReciprocalPath) -> Tuple{Vector{dtype(path)}, Vector{String}}

Get the position-label pairs of the ticks of a path.
"""
function ticks(path::ReciprocalPath)
    result = (dtype(path)[], String[])
    d = zero(dtype(path))
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
    title --> string(nameof(typeof(path)))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    coordinates = map(Tuple, path)
    @series begin
        seriestype := :scatter
        coordinates
    end
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
    coordinates
end

"""
    selectpath(brillouinzone::BrillouinZone, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; ends=nothing, atol::Real=atol, rtol::Real=rtol)
    selectpath(brillouinzone::BrillouinZone, points::NTuple{N, Number}...; labels=points, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N -> Tuple(ReciprocalPath, Vector{Int})
    selectpath(brillouinzone::BrillouinZone, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N -> Tuple(ReciprocalPath, Vector{Int})
    selectpath(brillouinzone::BrillouinZone, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N -> Tuple(ReciprocalPath, Vector{Int})
    selectpath(brillouinzone::BrillouinZone, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N -> Tuple(ReciprocalPath, Vector{Int})

Select a path from a `BrillouinZone`. Return a `ReciprocalPath` and the positions of the equivalent points in the `BrillouinZone`.

When `ends` is `nothing`, the start point will be included while the end point will be not for each segment except for the last whose both points will be included.
"""
@inline function selectpath(brillouinzone::BrillouinZone, contents::NamedTuple{(:points, :labels), <:Tuple{<:Tuple, <:Tuple}}; ends=nothing, atol::Real=atol, rtol::Real=rtol)
    return selectpath(brillouinzone, contents.points...; labels=contents.labels, ends=ends, atol=atol, rtol=rtol)
end
@inline function selectpath(brillouinzone::BrillouinZone, points::NTuple{N, Number}...; labels=points, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N
    return selectpath(brillouinzone, points; labels=labels, ends=ends, atol=atol, rtol=rtol)
end
@inline function selectpath(brillouinzone::BrillouinZone, points::Tuple{Vararg{NTuple{N, Number}}}; labels=points, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N
    return selectpath(brillouinzone, points2segments(points)...; labels=points2segments(Tuple(labels)), ends=ends, atol=atol, rtol=rtol)
end
@inline function selectpath(brillouinzone::BrillouinZone, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N
    return selectpath(brillouinzone, segments; labels=labels, ends=ends, atol=atol, rtol=rtol)
end
function selectpath(brillouinzone::BrillouinZone, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}}; labels=segments, ends=nothing, atol::Real=atol, rtol::Real=rtol) where N
    @assert length(brillouinzone.reciprocals)==N "selectpath error: mismatched number of reciprocals ($(length(brillouinzone.reciprocals)) v.s. $N)."
    isnothing(ends) || @assert fieldcount(typeof(segments))==length(ends) "selectpath error: mismatched number of segments."
    isnothing(ends) && (ends = ntuple(i->i<fieldcount(typeof(segments)) ? (true, segments[i].second==segments[i+1].first ? false : true) : (true, true), fieldcount(typeof(segments))))
    contents, indexes = ntuple(i->eltype(brillouinzone)[], fieldcount(typeof(segments))), Int[]
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
        positions, distances = Int[], dtype(brillouinzone)[]
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
    return ReciprocalPath{first(names(brillouinzone))}(contents, Tuple(labels)), indexes
end

const linemap = Dict(
    "Γ"=>(0//1,), "X"=>(1//2,), "Γ₁"=>(0//1,), "Γ₂"=>(1//1,), "X₁"=>(1//2,), "X₂"=>(-1//2,)
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
    "Γ"=>(0//1, 0//1), "X"=>(1//2, 0//1), "Y"=>(0//1, 1//2), "M"=>(1//2, 1//2),
    "X₁"=>(1//2, 0//1), "X₂"=>(-1//2, 0//1), "Y₁"=>(0//1, 1//2), "Y₂"=>(0//1, -1//2),
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
    "Γ"=>(0//1, 0//1), "K"=>(2//3, 1//3), "M"=>(1//2, 1//2),
    "K₁"=>(2//3, 1//3), "K₂"=>(1//3, 2//3), "K₃"=>(1//3, -1//3), "K₄"=>(-2//3, -1//3), "K₅"=>(-1//3, -2//3), "K₆"=>(-1//3, 1//3),
    "M₁"=>(1//2, 1//2), "M₂"=>(1//2, 0//1), "M₃"=>(0//1, -1//2), "M₄"=>(-1//2, -1//2), "M₅"=>(-1//2, 0//1), "M₆"=>(0//1, 1//2)
)
# High-symmetric point in the hexagonal Brillouin zone when the angle between the reciprocal vectors is 60°
const hexagon60°map = Dict(
    "Γ"=>(0//1, 0//1), "K"=>(1//3, 1//3), "M"=>(0//1, 1//2),
    "K₁"=>(1//3, 1//3), "K₂"=>(2//3, -1//3), "K₃"=>(1//3, -2//3), "K₄"=>(-1//3, -1//3), "K₅"=>(-2//3, 1//3), "K₆"=>(-1//3, 2//3),
    "M₁"=>(0//1, 1//2), "M₂"=>(1//2, 0//1), "M₃"=>(1//2, -1//2), "M₄"=>(0//1, -1//2), "M₅"=>(-1//2, 0//1), "M₆"=>(-1//2, 1//2)
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

# plot utilities
block = quote
    seriestype --> :path
    titlefontsize --> 10
    legend --> false
    minorgrid --> true
    xminorticks --> 10
    yminorticks --> 10
    xticks --> ticks(path)
    xlabel --> string(names(path)[1])
    [distance(path, i) for i=1:length(path)], data
end
@eval @recipe plot(path::ReciprocalPath, data::AbstractVector) = $block
@eval @recipe plot(path::ReciprocalPath, data::AbstractMatrix) = $block

@recipe function plot(path::ReciprocalPath, y::AbstractVector, data::AbstractMatrix)
    seriestype --> :heatmap
    titlefontsize --> 10
    xticks --> ticks(path)
    xlabel --> string(names(path)[1])
    xlims --> (0, distance(path))
    ylims --> (minimum(y), maximum(y))
    clims --> extrema(data)
    [distance(path, i) for i=1:length(path)], y, data
end

block = quote
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."
    x, y = xaxis(reciprocalspace), yaxis(reciprocalspace)
    Δx, Δy= x[2]-x[1], y[2]-y[1]
    seriestype --> :heatmap
    titlefontsize --> 10
    aspect_ratio --> :equal
    xlims --> (x[1]-Δx, x[end]+Δx)
    ylims --> (y[1]-Δy, y[end]+Δy)
    clims --> extrema(data)
    xlabel --> string(names(reciprocalspace)[1], "₁")
    ylabel --> string(names(reciprocalspace)[1], "₂")
    x, y, data
end
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractMatrix) = $block
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractMatrix) = $block

setup(expr::Expr) = quote
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
    seriestype := :heatmap
    xlims := (minimum(clims), maximum(clims))
    xlabel := ""
    ylims := (0, 1)
    yticks := (0:1, ("", ""))
    ylabel := ""
    LinRange(clims..., 100), [0, 1], [LinRange(clims..., 100)'; LinRange(clims..., 100)']
end
@eval @recipe plot(path::ReciprocalPath, y::AbstractVector, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(path, y, data[:, :, i])))
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))

# save utilities
function save(filename::AbstractString, path::ReciprocalPath, data::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}})
    @assert length(path)==size(data)[1] "save error: mismatched size of path and data."
    open(filename, "w") do f
        writedlm(f, [matrix(path) data])
    end
end
function save(filename::AbstractString, path::ReciprocalPath, y::AbstractVector, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}})
    @assert size(data)[1:2]==(length(y), length(path)) "save error: mismatched size of path, y and data."
    open(filename, "w") do f
        x = matrix(kron(path, ones(length(y))))
        y = kron(ones(length(path)), y)
        z = reshape(data, prod(size(data)[1:2]), :)
        writedlm(f, [x y z])
    end
end
function save(filename::AbstractString, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}})
    @assert length(reciprocalspace)==prod(size(data)[1:2]) "save error: mismatched size of reciprocal space and data."
    open(filename, "w") do f
        writedlm(f, [matrix(reciprocalspace) reshape(data, prod(size(data)[1:2]), :)])
    end
end
function matrix(vs::AbstractVector{<:AbstractVector})
    result = zeros(eltype(eltype(vs)), length(vs), length(vs[1]))
    for (i, v) in enumerate(vs)
        result[i, :] = v
    end
    return result
end

end #module
