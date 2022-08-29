module Spatials

using Base.Iterators: flatten, product
using LinearAlgebra: cross, dot, norm
using NearestNeighbors: KDTree, inrange, knn
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series
using StaticArrays: SVector

using ..QuantumNumbers: AbelianNumbers, Momentum, Momentum₁, Momentum₂, Momentum₃, periods
using ...Prerequisites: atol, rtol, Float
using ...Prerequisites: CompositeDict
using ...Prerequisites.Traits: efficientoperations
using ...Prerequisites.VectorSpaces: VectorSpace, SimpleNamedVectorSpace, VectorSpaceStyle, VectorSpaceCartesian

import ..Essentials: dtype, kind, reset!
import ...Interfaces: decompose, dimension
import ...Prerequisites.Traits: contentnames, getcontent
import ...Prerequisites.VectorSpaces: shape

export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, Segment, Translations, bonds!, bonds, icoordinate, isintracell, nneighbor, rcoordinate, @translations_str
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

"""
    isintratriangle(p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
        vertexes::NTuple{3, Bool}=(true, true, true), edges::NTuple{3, Bool}=(true, true, true), atol::Real=atol, rtol::Real=rtol
        ) -> Bool

Judge whether a point belongs to the interior of a triangle whose vertexes are `p₁`, 'p₂' and `p₃` with the give tolerance. `vertexes` and `edges` define whether the interior should contain the vertexes or edges, respectively.
!!! note
    1. The vertexes are in the order (p₁, p₂, p₃) and the edges are in the order (p1p2, p2p3, p3p1).
    2. The edges do not contain the vertexes.
"""
function isintratriangle(p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number}, p₃::AbstractVector{<:Number};
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
    isonline(p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number};
        ends::Tuple{Bool, Bool}=(true, true), atol::Real=atol, rtol::Real=rtol
        ) -> Bool

Judge whether a point is on a line segment whose end points are `p₁` and `p₂` with the given tolerance. `ends` defines whether the line segment should contain its ends.
"""
function isonline(p::AbstractVector{<:Number}, p₁::AbstractVector{<:Number}, p₂::AbstractVector{<:Number};
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
    reciprocals(vectors::AbstractVector{AbstractVector{<:Number}}) -> Vector{Vector{<:Number}}

Get the reciprocals dual to the input vectors.
"""
function reciprocals(vectors::AbstractVector{<:AbstractVector{<:Number}})
    @assert length(vectors)<4 "reciprocals error: the number of input vectors should not be greater than 3."
    @assert mapreduce(v->length(v)∈(1, 2, 3), &, vectors, init=true) "reciprocals error: all input vectors must be 1, 2 or 3 dimensional."
    datatype = promote_type(Float, eltype(eltype(vectors)))
    result = Vector{datatype}[]
    if length(vectors) == 1
        push!(result, 2*convert(datatype, pi)/mapreduce(vi->vi^2, +, vectors[1])*vectors[1])
    elseif length(vectors) == 2
        v₁, v₂ = vectors[1], vectors[2]
        @assert length(v₁)==length(v₂) "reciprocals error: mismatched length of input vectors."
        if length(v₁) == 2
            det = 2*convert(datatype, pi) / (v₁[1]*v₂[2]-v₁[2]*v₂[1])
            push!(result, [det*v₂[2], -det*v₂[1]])
            push!(result, [-det*v₁[2], det*v₁[1]])
        else
            v₃ = SVector{3, datatype}(v₁[2]*v₂[3]-v₁[3]*v₂[2], v₁[3]*v₂[1]-v₁[1]*v₂[3], v₁[1]*v₂[2]-v₁[2]*v₂[1])
            V = 2*convert(datatype, pi) / volume(v₁, v₂, v₃)
            push!(result, [v₂[2]*v₃[3]*V-v₂[3]*v₃[2]*V, v₂[3]*v₃[1]*V-v₂[1]*v₃[3]*V, v₂[1]*v₃[2]*V-v₂[2]*v₃[1]*V])
            push!(result, [v₃[2]*v₁[3]*V-v₃[3]*v₁[2]*V, v₃[3]*v₁[1]*V-v₃[1]*v₁[3]*V, v₃[1]*v₁[2]*V-v₃[2]*v₁[1]*V])
        end
    elseif length(vectors) == 3
        v₁, v₂, v₃ = vectors
        V = 2*convert(datatype, pi) / volume(v₁, v₂, v₃)
        push!(result, [v₂[2]*v₃[3]*V-v₂[3]*v₃[2]*V, v₂[3]*v₃[1]*V-v₂[1]*v₃[3]*V, v₂[1]*v₃[2]*V-v₂[2]*v₃[1]*V])
        push!(result, [v₃[2]*v₁[3]*V-v₃[3]*v₁[2]*V, v₃[3]*v₁[1]*V-v₃[1]*v₁[3]*V, v₃[1]*v₁[2]*V-v₃[2]*v₁[1]*V])
        push!(result, [v₁[2]*v₂[3]*V-v₁[3]*v₂[2]*V, v₁[3]*v₂[1]*V-v₁[1]*v₂[3]*V, v₁[1]*v₂[2]*V-v₁[2]*v₂[1]*V])
    end
    return result
end

"""
    rotate(cluster::AbstractMatrix{<:Number}, angle::Number;
        axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))
        ) -> Matrix{<:Number}

Get a rotated cluster of the original one by a certain angle around an axis.

The axis is determined by a point it gets through (`nothing` can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.
!!! note
    1. The result is given by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
    2. Only 2 and 3 dimensional vectors can be rotated.
    3. When the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.
"""
function rotate(cluster::AbstractMatrix{<:Number}, angle::Number;
        axis::Tuple{Union{AbstractVector{<:Number}, Nothing}, Tuple{<:Number, <:Number}}=(nothing, (0, 0))
        )
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
    minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Int=1; coordination::Int=8) -> Vector{Float}

Use kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.

When the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds across the periodic boundaries, the cluster will be pre-translated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The `coordination` parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain `Inf`. This is a sign that you may need a larger `coordination`. Another situation that `Inf` appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.
"""
function minimumlengths(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, nneighbor::Int=1; coordination::Int=8)
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
    return result
end

"""
    Neighbors{K, V<:Number} <: CompositeDict{K, V}

Neighbor vs. bond length maps.
"""
struct Neighbors{K, V<:Number} <: CompositeDict{K, V}
    contents::Dict{K, V}
end
@inline Neighbors(pairs...) = Neighbors(Dict(pairs...))
@inline Neighbors(lengths::Vector{<:Number}) = Neighbors(i=>length for (i, length) in enumerate(lengths))
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
    Translations{N} <: VectorSpace{NTuple{N, Int}}

A set of translations. The boundary conditions along each translation direction are also included.
"""
struct Translations{N} <: VectorSpace{NTuple{N, Int}}
    ranges::NTuple{N, Int}
    boundaries::NTuple{N, Char}
    function Translations(ranges::NTuple{N, Int}, boundaries::NTuple{N, Char}) where N
        @assert all(map(>(0), ranges)) "Translations error: ranges must be positive."
        @assert all(map(in(('P', 'O', 'p', 'o')), boundaries)) "Translations error: boundary conditions must be either 'P'/'p' for 'periodic' or 'O'/'o' for 'open'."
        new{N}(ranges, map(uppercase, boundaries))
    end
end
@inline VectorSpaceStyle(::Type{<:Translations}) = VectorSpaceCartesian()
@inline shape(translations::Translations) = map(i->0:i-1, translations.ranges)
@inline Base.CartesianIndex(translation::NTuple{N, Int}, ::Translations{N}) where N = CartesianIndex(translation)
@inline Tuple(index::CartesianIndex{N}, ::Translations{N}) where N = Tuple(index)
function Base.show(io::IO, translations::Translations{N}) where N
    for i = 1:N
        i>1 && @printf io "%s" '-'
        @printf io "%s%s" translations.ranges[i] translations.boundaries[i]
    end
end

"""
    translations"n₁P/O-n₁P/O-..." -> Translations

Construct a set of translations from a string literal.
"""
macro translations_str(str::String)
    patterns = split(str, "-")
    ranges, boundaries = Int[], Char[]
    for pattern in patterns
        push!(ranges, Meta.parse(pattern[1:end-1]))
        push!(boundaries, pattern[end])
    end
    return Translations(Tuple(ranges), Tuple(boundaries))
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
@inline Base.show(io::IO, bond::Bond) = @printf io "Bond(%s, %s)" bond.kind join(map(string, bond.points), ", ")

"""
    Bond(point::Point)
    Bond(kind::Integer, point::Point, points::Point...)

Construct a bond.
"""
@inline Bond(point::Point) = Bond(0, point)
@inline Bond(kind::Integer, point::Point, points::Point...) = Bond(kind, [point, points...])

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
    return bond[1].rcoordinate-bond[2].rcoordinate
end

"""
    icoordinate(bond::Bond) -> SVector

Get the icoordinate of the bond.
"""
@inline function icoordinate(bond::Bond)
    @assert length(bond)∈(1, 2) "icoordinate error: not supported for a generic bond which contains $(length(bond)) points."
    length(bond)==1 && return bond[1].icoordinate
    return bond[1].icoordinate-bond[2].icoordinate
end

"""
    isintracell(bond::Bond) -> Bool

Judge whether a bond is intra the unit cell of a lattice.
"""
@inline isintracell(bond::Bond) = all(point->isapprox(norm(point.icoordinate), 0.0, atol=atol, rtol=rtol), bond)

"""
    AbstractLattice{N, D<:Number}

Abstract type of a unitcell-described lattice.

It should have the following contents:
- `name::Symbol`: the name of the lattice
- `coordinates::Matrix{D}`: the coordinates of the lattice
- `vectors::Vector{SVector{N, D}}`: the translation vectors of the lattice
- `reciprocals::Vector{SVector{N, D}}`: the reciprocals of the lattice
"""
abstract type AbstractLattice{N, D<:Number} end
@inline contentnames(::Type{<:AbstractLattice}) = (:name, :coordinates, :vectors, :reciprocals)
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
    Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Int=8)

Get the neighbor vs. bond length map of a lattice up to the `nneighbor`th order.
"""
@inline Neighbors(lattice::AbstractLattice, nneighbor::Integer; coordination::Int=8) = Neighbors(minimumlengths(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), nneighbor; coordination=coordination))

"""
    bonds(lattice::AbstractLattice, nneighbor::Int; coordination::Int=8) -> Vector{Bond{Int, Point{dimension(lattice), dtype(lattice)}}}
    bonds(lattice::AbstractLattice, neighbors::Neighbors) -> Vector{Bond{keytype(neighbors), Point{dimension(lattice), dtype(lattice)}}}

Get the required bonds of a lattice.
"""
@inline bonds(lattice::AbstractLattice, nneighbor::Int; coordination::Int=8) = bonds!(Bond{Int, Point{dimension(lattice), dtype(lattice)}}[], lattice, nneighbor; coordination=coordination)
@inline bonds(lattice::AbstractLattice, neighbors::Neighbors) = bonds!(Bond{keytype(neighbors), Point{dimension(lattice), dtype(lattice)}}[], lattice, neighbors)

"""
    bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Int; coordination::Int=8)
    bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors) -> typeof(bonds)

Get the required bonds of a lattice and append them to the input bonds.
"""
@inline bonds!(bonds::Vector, lattice::AbstractLattice, nneighbor::Int; coordination::Int=8) = bonds!(bonds, lattice, Neighbors(lattice, nneighbor; coordination=coordination))
function bonds!(bonds::Vector, lattice::AbstractLattice, neighbors::Neighbors)
    origin = zero(SVector{dimension(lattice), dtype(lattice)})
    for (i, coordinate) in enumerate(lattice)
        push!(bonds, Bond(0, Point(i, coordinate, origin)))
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
    @recipe plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true)

Define the recipe for the visualization of a lattice.
"""
@recipe function plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true)
    title := String(getcontent(lattice, :name))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    @series begin
        seriestype := :scatter
        coordinates = NTuple{dimension(lattice), dtype(lattice)}[]
        for i = 1:length(lattice)
            bond = Bond(0, Point(i, lattice[i], zero(lattice[i])))
            filter(bond) && push!(coordinates, Tuple(lattice[i]))
        end
        coordinates
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
    Lattice{N, D<:Number} <: AbstractLattice{N, D}

Simplest lattice.

A simplest lattice can be constructed from its coordinates and translation vectors.
"""
struct Lattice{N, D<:Number} <: AbstractLattice{N, D}
    name::Symbol
    coordinates::Matrix{D}
    vectors::Vector{SVector{N, D}}
    reciprocals::Vector{SVector{N, D}}
    function Lattice{N}(name::Symbol, coordinates::Matrix{<:Number}, vectors::Vector{<:AbstractVector{<:Number}}) where N
        @assert N==size(coordinates, 1) "Lattice error: shape mismatched."
        datatype = promote_type(Float, eltype(coordinates), eltype(eltype(vectors)))
        coordinates = convert(Matrix{datatype}, coordinates)
        vectors = convert(Vector{SVector{N, datatype}}, vectors)
        recipls = convert(Vector{SVector{N, datatype}}, reciprocals(vectors))
        new{N, datatype}(name, coordinates, vectors, recipls)
    end
end

"""
    Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{Vector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    Lattice(coordinates::Vector{<:Number}...; name::Symbol=:lattice, vectors::Union{Vector{<:AbstractVector{<:Number}}, Nothing}=nothing)

Construct a lattice.
"""
function Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{Vector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    isnothing(vectors) && (vectors = SVector{0, SVector{N, eltype(eltype(coordinates))}}())
    coordinates = [coordinates[j][i] for i=1:N, j=1:length(coordinates)]
    return Lattice{N}(name, coordinates, vectors)
end
function Lattice(coordinates::Vector{<:Number}...; name::Symbol=:lattice, vectors::Union{Vector{<:AbstractVector{<:Number}}, Nothing}=nothing)
    coordinates = hcat(coordinates...)
    isnothing(vectors) && (vectors = SVector{size(coordinates)[1], eltype(coordinates)}[])
    return Lattice{size(coordinates)[1]}(name, coordinates, vectors)
end

"""
    Lattice(lattice::Lattice, translations::Translations)

Construct a lattice from the translations of another.
"""
function Lattice(lattice::Lattice, translations::Translations)
    name = Symbol(@sprintf "%s(%s)" lattice.name translations)
    coordinates = tile(lattice.coordinates, lattice.vectors, translations)
    vectors = SVector{dimension(lattice), dtype(lattice)}[]
    for (i, vector) in enumerate(lattice.vectors)
        translations.boundaries[i]=='P' && push!(vectors, vector*translations.ranges[i])
    end
    return Lattice{dimension(lattice)}(name, coordinates, vectors)
end

"""
    Segment{S} <: AbstractVector{S}

A segment.
"""
struct Segment{S} <: AbstractVector{S}
    start::S
    stop::S
    length::Int
    ends::Tuple{Bool, Bool}
end
@inline Base.:(==)(s₁::Segment, s₂::Segment) = ==(efficientoperations, s₁, s₂)
@inline Base.isequal(s₁::Segment, s₂::Segment) = isequal(efficientoperations, s₁, s₂)
@inline Base.size(segment::Segment) = (segment.length,)
function Base.getindex(segment::Segment, i::Integer)
    length = segment.length+count(isequal(false), segment.ends)-1
    step = (segment.stop-segment.start)/length
    start = segment.ends[1] ? segment.start : segment.start+step
    return start+(i-1)*step
end
function iterate(segment::Segment)
    segment.length==0 && return
    length = segment.length+count(isequal(false), segment.ends)-1
    step = (segment.stop-segment.start)/length
    start = segment.ends[1] ? segment.start : segment.start+step
    return start, (1, start, step)
end
function iterate(segment::Segment, state)
    i, middle, step = state
    i==segment.length && return
    middle = middle+step
    return middle, (i+1, middle, step)
end
function Base.show(io::IO, segment::Segment{<:Number})
    left, right = (segment.ends[1] ? "[" : "("), (segment.ends[2] ? "]" : ")")
    @printf io "%s%s, %s%s" left segment.start segment.stop right
end
function Base.show(io::IO, segment::Segment)
    left, right = (segment.ends[1] ? "[" : "("), (segment.ends[2] ? "]" : ")")
    @printf io "%sp₁, p₂%s with p₁=%s and p₂=%s" left right segment.start segment.stop
end

"""
    Segment(start::Number, stop::Number, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    Segment(start::AbstractVector, stop::AbstractVector, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N

Construct a segment.
"""
function Segment(start::Number, stop::Number, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = promote_type(typeof(start), typeof(stop), Float64)
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::AbstractVector, stop::AbstractVector, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    @assert Base.length(start)==Base.length(stop) "Segment error: start and stop should have equal length."
    dtype = SVector{Base.length(start), promote_type(eltype(start), eltype(stop), Float64)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = SVector{N, promote_type(eltype(start), eltype(stop), Float64)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end

"""
    ReciprocalSpace{K, P} <: SimpleNamedVectorSpace{K, P}

Abstract type of reciprocal spaces.
"""
abstract type ReciprocalSpace{K, P} <: SimpleNamedVectorSpace{K, P} end

@inline vectorconvert(reciprocals::Vector{<:SVector}) = reciprocals
@inline vectorconvert(reciprocals::AbstractVector) = convert(Vector{SVector{length(first(reciprocals)), eltype(eltype(reciprocals))}}, reciprocals)
@inline vectorconvert(reciprocals::AbstractVector{<:SVector}) = convert(Vector{SVector{length(eltype(reciprocals)), eltype(eltype(reciprocals))}}, reciprocals)
"""
    BrillouinZone{K, P<:Momentum, S<:SVector} <: ReciprocalSpace{K, P}

The Brillouin zone of a lattice.
"""
struct BrillouinZone{K, P<:Momentum, S<:SVector} <: ReciprocalSpace{K, P}
    reciprocals::Vector{S}
    momenta::AbelianNumbers{P}
    function BrillouinZone{K}(reciprocals::AbstractVector, momenta::AbelianNumbers{<:Momentum}) where K
        @assert isa(K, Symbol) "BrillouinZone error: K must be a Symbol."
        reciprocals = vectorconvert(reciprocals)
        new{K, eltype(momenta), eltype(reciprocals)}(reciprocals, momenta)
    end
end
@inline contentnames(::Type{<:BrillouinZone}) = (:reciprocals, :content)
@inline getcontent(bz::BrillouinZone, ::Val{:content}) = bz.momenta

"""
    BrillouinZone(reciprocals::AbstractVector, nk::Int)
    BrillouinZone(::Type{P}, reciprocals::AbstractVector) where {P<:Momentum}
    BrillouinZone(reciprocals::AbstractVector, momenta::AbelianNumbers{<:Momentum})

    BrillouinZone{K}(reciprocals::AbstractVector, nk::Int) where K
    BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector) where {K, P<:Momentum}
    BrillouinZone{K}(reciprocals::AbstractVector, momenta::AbelianNumbers{<:Momentum}) where K

Construct a Brillouin zone.
"""
@inline BrillouinZone(reciprocals::AbstractVector, nk::Int) = BrillouinZone{:k}(reciprocals, nk)
@inline BrillouinZone(::Type{P}, reciprocals::AbstractVector) where {P<:Momentum} = BrillouinZone{:k}(P, reciprocals)
@inline BrillouinZone(reciprocals::AbstractVector, momenta::AbelianNumbers{<:Momentum}) = BrillouinZone{:k}(reciprocals, momenta)
@inline function BrillouinZone{K}(reciprocals::AbstractVector, nk::Int) where K
    length(reciprocals)==1 && return BrillouinZone{K}(Momentum₁{nk}, reciprocals)
    length(reciprocals)==2 && return BrillouinZone{K}(Momentum₂{nk, nk}, reciprocals)
    length(reciprocals)==3 && return BrillouinZone{K}(Momentum₃{nk, nk, nk}, reciprocals)
    error("BrillouinZone error: only 1d, 2d and 3d are supported.")
end
@inline function BrillouinZone{K}(::Type{P}, reciprocals::AbstractVector) where {K, P<:Momentum}
    contents = Vector{P}(undef, prod(periods(P)))
    for (i, index) in enumerate(CartesianIndices(reverse(periods(P))))
        contents[i] = P(reverse((index-oneunit(index)).I)...)
    end
    momenta = AbelianNumbers('C', contents)
    return BrillouinZone{K}(reciprocals, momenta)
end

"""
    dtype(bz::BrillouinZone)
    dtype(::Type{<:BrillouinZone{K, <:Momentum, S} where K}) where {S<:SVector}

Get the data type of a Brillouin zone.
"""
@inline dtype(bz::BrillouinZone) = dtype(typeof(bz))
@inline dtype(::Type{<:BrillouinZone{K, <:Momentum, S} where K}) where {S<:SVector} = eltype(S)

"""
    ReciprocalZone{K, S<:SVector, V<:Number} <: ReciprocalSpace{K, S}

A zone in the reciprocal space.
"""
struct ReciprocalZone{K, S<:SVector, V<:Number} <: ReciprocalSpace{K, S}
    momenta::Vector{S}
    reciprocals::Vector{S}
    bounds::Vector{Segment{V}}
    volume::V
    function ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Vector{<:Segment}) where K
        @assert isa(K, Symbol) "ReciprocalZone error: K must be a Symbol."
        @assert length(reciprocals)==length(bounds) "ReciprocalZone error: mismatched number of reciprocals and its boundaries."
        reciprocals = vectorconvert(reciprocals)
        momenta = zeros(SVector{length(eltype(reciprocals)), promote_type(eltype(eltype(reciprocals)), eltype(eltype(bounds)))}, prod(map(length, bounds)))
        for (i, index) in enumerate(product(reverse(bounds)...))
            index = reverse(index)
            for j = 1:length(index)
                momenta[i] += reciprocals[j]*index[j]
            end
        end
        ratio = one(eltype(eltype(momenta)))
        for i = 1:length(bounds)
            ratio = ratio*(bounds[i].stop-bounds[i].start)
        end
        v = ratio*volume(reciprocals)
        return new{K, eltype(momenta), typeof(v)}(momenta, reciprocals, bounds, v)
    end
end
@inline contentnames(::Type{<:ReciprocalZone}) = (:content, :reciprocals, :bounds, :volume)
@inline getcontent(rz::ReciprocalZone, ::Val{:content}) = rz.momenta

"""
    ReciprocalZone(reciprocals::AbstractVector; length::Int=100)
    ReciprocalZone(reciprocals::AbstractVector, bounds::Segment...)
    ReciprocalZone(reciprocals::AbstractVector, bounds::Union{Tuple{Vararg{Segment}}, Vector{<:Segment}})

    ReciprocalZone{K}(reciprocals::AbstractVector; length::Int=100) where K
    ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Segment...) where K
    ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Union{Tuple{Vararg{Segment}}, Vector{<:Segment}}) where K

Construct a rectangular zone in the reciprocal space.
"""
@inline ReciprocalZone(reciprocals::AbstractVector; length::Int=100) = ReciprocalZone{:k}(vectorconvert(reciprocals); length=length)
@inline ReciprocalZone(reciprocals::AbstractVector, bounds::Segment...) = ReciprocalZone{:k}(vectorconvert(reciprocals), bounds)
@inline ReciprocalZone(reciprocals::AbstractVector, bounds::Union{Tuple{Vararg{Segment}}, Vector{<:Segment}}) = ReciprocalZone{:k}(vectorconvert(reciprocals), bounds)
@inline function ReciprocalZone{K}(reciprocals::AbstractVector; length::Int=100) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), ntuple(i->Segment(-1//2, 1//2, length, ends=(true, false)), Base.length(reciprocals)))
end
@inline ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Segment...) where K = ReciprocalZone{K}(vectorconvert(reciprocals), bounds)
@inline ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Tuple{Vararg{Segment}}) where K = ReciprocalZone{K}(vectorconvert(reciprocals), collect(bounds))

"""
    ReciprocalZone(bz::BrillouinZone)

Construct a reciprocal zone from a Brillouin zone.
"""
@inline function ReciprocalZone(bz::BrillouinZone{K, P}) where {K, P<:Momentum}
    return ReciprocalZone{K}(bz.reciprocals, map(length->Segment(0, 1, length, ends=(true, false)), periods(P)))
end

"""
    ReciprocalPath{K, S<:SVector} <: ReciprocalSpace{K, S}

A path in the reciprocal space.
"""
struct ReciprocalPath{K, S<:SVector} <: ReciprocalSpace{K, S}
    momenta::Vector{S}
    function ReciprocalPath{K}(momenta::AbstractVector) where K
        @assert isa(K, Symbol) "ReciprocalPath error: K must be a Symbol."
        momenta = vectorconvert(momenta)
        new{K, eltype(momenta)}(momenta)
    end
end
@inline contentnames(::Type{<:ReciprocalPath}) = (:content,)
@inline getcontent(rp::ReciprocalPath, ::Val{:content}) = rp.momenta

"""
    ReciprocalPath(momenta::AbstractVector)
    ReciprocalPath(reciprocals::AbstractVector, segments::Segment{<:SVector}...)
    ReciprocalPath(reciprocals::AbstractVector, segments::Tuple{Vararg{Segment{<:SVector}}})

    ReciprocalPath{K}(momenta::AbstractVector) where K
    ReciprocalPath{K}(reciprocals::AbstractVector, segments::Segment{<:SVector}...) where K
    ReciprocalPath{K}(reciprocals::AbstractVector, segments::Tuple{Vararg{Segment{<:SVector}}}) where K

Construct a path in the reciprocal space.
"""
@inline ReciprocalPath(momenta::AbstractVector) = ReciprocalPath{:k}(momenta)
@inline ReciprocalPath(reciprocals::AbstractVector, segments::Segment{<:SVector}...) = ReciprocalPath{:k}(vectorconvert(reciprocals), segments)
@inline ReciprocalPath(reciprocals::AbstractVector, segments::Tuple{Vararg{Segment{<:SVector}}}) = ReciprocalPath{:k}(vectorconvert(reciprocals), segments)
@inline ReciprocalPath{K}(reciprocals::AbstractVector, segments::Segment{<:SVector}...) where K = ReciprocalPath{K}(vectorconvert(reciprocals), segments)
@inline function ReciprocalPath{K}(reciprocals::AbstractVector, segments::Tuple{Vararg{Segment{<:SVector}}}) where K
    @assert length(reciprocals)==length(eltype(eltype(segments))) "ReciprocalPath error: mismatched number of reciprocals."
    reciprocals = vectorconvert(reciprocals)
    momenta = zeros(SVector{length(eltype(reciprocals)), promote_type(eltype(eltype(reciprocals)), eltype(eltype(eltype(segments))))}, sum(map(length, segments)))
    count = 1
    for segment in segments
        for index in segment
            for i = 1:length(index)
                momenta[count] += reciprocals[i]*index[i]
            end
            count += 1
        end
    end
    return ReciprocalPath{K}(momenta)
end

"""
    ReciprocalPath(reciprocals::AbstractVector, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...;
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {N, M}
    ReciprocalPath(reciprocals::AbstractVector, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}};
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {N, M}

    ReciprocalPath{K}(reciprocals::AbstractVector, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...;
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {K, N, M}
    ReciprocalPath{K}(reciprocals::AbstractVector, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}};
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {K, N, M}

Construct a path in the reciprocal space.
"""
@inline function ReciprocalPath(reciprocals::AbstractVector, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...;
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {N, M}
    return ReciprocalPath{:k}(reciprocals, segments, length=length, ends=ends)
end
function ReciprocalPath(reciprocals::AbstractVector, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}};
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {N, M}
    return ReciprocalPath{:k}(reciprocals, segments, length=length, ends=ends)
end
@inline function ReciprocalPath{K}(reciprocals::AbstractVector, segments::Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}...;
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {K, N, M}
    return ReciprocalPath{K}(reciprocals, segments, length=length, ends=ends)
end
function ReciprocalPath{K}(reciprocals::AbstractVector, segments::Tuple{Vararg{Pair{<:NTuple{N, Number}, <:NTuple{N, Number}}}};
        length::Union{Int, NTuple{M, Int}}=100,
        ends::Union{NTuple{2, Bool}, NTuple{M, NTuple{2, Bool}}}=(true, false)
        ) where {K, N, M}
    @assert Base.length(reciprocals)==N "ReciprocalPath error: mismatched number of reciprocals ($(Base.length(reciprocals)) v.s. $N)."
    (isa(length, Int) && isa(ends, NTuple{2, Bool})) || @assert fieldcount(typeof(segments))==M "ReciprocalPath error: mismatched number of segments."
    isa(length, Int) && (length = ntuple(i->length, Val(fieldcount(typeof(segments)))))
    isa(ends, NTuple{2, Bool}) && (ends = ntuple(i->ends, Val(fieldcount(typeof(segments)))))
    segments = ntuple(i->Segment(segments[i].first, segments[i].second, length[i]; ends=ends[i]), Val(fieldcount(typeof(segments))))
    return ReciprocalPath{K}(reciprocals, segments)
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
    return ntuple(i->linemap[points[i]]=>linemap[points[i+1]], length(points)-1)
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
    return ntuple(i->rectanglemap[points[i]]=>rectanglemap[points[i+1]], length(points)-1)
end

const hexagon120°map = Dict(
    "Γ"=>(0//1, 0//1), "K"=>(2//3, 1//3), "M"=>(1//2, 1//2),
    "K₁"=>(2//3, 1//3), "K₂"=>(1//3, 2//3), "K₃"=>(1//3, -1//3), "K₄"=>(-2//3, -1//3), "K₅"=>(-1//3, -2//3), "K₆"=>(-1//3, 1//3),
    "M₁"=>(1//2, 1//2), "M₂"=>(1//2, 0//1), "M₃"=>(0//1, -1//2), "M₄"=>(-1//2, -1//2), "M₅"=>(-1//2, 0//1), "M₆"=>(0//1, 1//2)
)
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
    return ntuple(i->map[points[i]]=>map[points[i+1]], length(points)-1)
end

end #module
