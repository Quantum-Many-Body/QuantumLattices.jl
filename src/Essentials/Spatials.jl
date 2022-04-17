module Spatials

using LinearAlgebra: norm, dot, cross
using StaticArrays: SVector
using Printf: @printf, @sprintf
using NearestNeighbors: KDTree, knn, inrange
using Base.Iterators: flatten, product
using RecipesBase: RecipesBase, @recipe, @series
using ..QuantumOperators: OperatorUnit
using ..QuantumNumbers: AbelianNumbers, Momentum, Momentum₁, Momentum₂, Momentum₃, periods
using ...Prerequisites: atol, rtol, Float
using ...Prerequisites.Combinatorics: Combinations
using ...Prerequisites.Traits: efficientoperations, getcontent
using ...Prerequisites.SimpleTrees: SimpleTree, simpletreedepth, isleaf
using ...Prerequisites.VectorSpaces: VectorSpace, SimpleNamedVectorSpace, VectorSpaceStyle, VectorSpaceCartesian

import ..Essentials: dtype, kind, reset!
import ...Interfaces: decompose, rank, dimension, expand
import ...Prerequisites.Traits: contentnames, getcontent
import ...Prerequisites.VectorSpaces: shape

export distance, azimuthd, azimuth, polard, polar, volume, isparallel, isonline, isintratriangle, issubordinate
export reciprocals, translate, rotate, tile, minimumlengths, intralinks, interlinks, Translations, @translations_str
export AbstractPID, PID, CPID, AbstractBond, Point, Bond, pidtype, rcoord, icoord, isintracell
export AbstractLattice, Lattice, SuperLattice, Cylinder, LatticeIndex, LatticeBonds, Bonds
export nneighbor, bonds!, bonds, latticetype, bondtypes, latticebondsstructure
export allbonds, zerothbonds, insidebonds, acrossbonds, intrabonds, interbonds
export Segment, ReciprocalSpace, BrillouinZone, ReciprocalZone, ReciprocalPath
export linemap, rectanglemap, hexagon120°map, hexagon60°map, @line_str, @rectangle_str, @hexagon_str

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
    polard(v::AbstractVector{<:Number}) -> Number

Get the polar angle in degrees of a vector.
"""
function polard(v::AbstractVector{<:Number})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acosd(v[3]/norm(v))
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
    issubordinate(rcoord::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}; atol::Real=atol, rtol::Real=rtol) -> Bool

Judge whether a coordinate belongs to a lattice defined by `vectors` with the given tolerance.
"""
function issubordinate(rcoord::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}; atol::Real=atol, rtol::Real=rtol)
    @assert length(rcoord)∈(1, 2, 3) "issubordinate error: only 1, 2 and 3 dimensional coordinates are supported."
    @assert length(vectors)∈(1, 2, 3) "issubordinate error: the number of input basis vectors must be 1, 2 or 3."
    fapprox = xi->isapprox(round(xi), xi, atol=atol, rtol=rtol)
    if length(vectors) == 1
        result = mapreduce(fapprox, &, decompose(rcoord, vectors[1]))
    elseif length(vectors) == 2
        result = mapreduce(fapprox, &, decompose(rcoord, vectors[1], vectors[2]))
    else
        result = mapreduce(fapprox, &, decompose(rcoord, vectors[1], vectors[2], vectors[3]))
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
    translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number}) -> Matrix{vector|>eltype}

Get the translated cluster of the original one by a vector.
"""
@inline translate(cluster::AbstractMatrix{<:Number}, vector::AbstractVector{<:Number}) = cluster .+ reshape(vector, (vector|>length, 1))

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
    intralinks(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, neighbors::Dict{Int, <:Number},
        maxtranslations::NTuple{N, Int}=ntuple(i->length(neighbors), length(vectors))
        ) where N -> Vector{Tuple{Int, Int, Int, Vector{<:Number}}}

Use kdtree to get the intracluster nearest neighbors.

As is similar to [`minimumlengths`](@ref), when `vectors` is nonempty, the cluster assumes periodic boundaries. `neighbors` provides the map between the bond length and the order of nearest neighbors. Note only those with the lengths present in `neighbors` will be included in the result. `maxtranslations` determines the maximum number of translations along those directions specified by `vectors` when the tiled supercluster is constructed (See [`minimumlengths`](@ref) for the explanation of the method for periodic lattices).
"""
function intralinks(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, neighbors::Dict{Int, <:Number},
        maxtranslations::NTuple{N, Int}=ntuple(i->length(neighbors), length(vectors))
        ) where N
    @assert length(vectors)==N "intralinks error: mismatched shape of input vectors and maxtranslations."
    datatype = promote_type(eltype(cluster), eltype(eltype(vectors)))
    result = Tuple{Int, Int, Int, Vector{datatype}}[]
    length(neighbors)==0 && return result
    translations = reshape(product((-nnb:nnb for nnb in maxtranslations)...)|>collect, :)
    for translation in translations
        if length(translation)>0 && mapreduce(≠(0), |, translation)
            remove = ntuple(i->-translation[i], length(translation))
            filter!(≠(remove), translations)
        end
    end
    sort!(translations, by=norm)
    supercluster = tile(cluster, vectors, translations)
    disps = tile(zero(cluster), vectors, translations)
    for (i, indices) in enumerate(inrange(KDTree(convert(Matrix{Float}, supercluster)), convert(Matrix{Float}, cluster), max(values(neighbors)...)+atol, true))
        for j in indices
            if i < j
                dist = zero(datatype)
                for k = 1:size(cluster, 1)
                    dist = dist + (supercluster[k, j]-cluster[k, i])^2
                end
                dist = sqrt(dist)
                for (nb, len) in neighbors
                    if isapprox(len, dist, atol=atol, rtol=rtol)
                        push!(result, (nb, i, (j-1)%size(cluster, 2)+1, disps[:, j]))
                        break
                    end
                end
            end
        end
    end
    return result
end

"""
    interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Dict{Int, <:Number}) -> Vector{Tuple{Int, Int, Int}}

Use kdtree to get the intercluster nearest neighbors.
"""
function interlinks(cluster₁::AbstractMatrix{<:Number}, cluster₂::AbstractMatrix{<:Number}, neighbors::Dict{Int, <:Number})
    @assert size(cluster₁, 1)==size(cluster₂, 1) "interlinks error: mismatched space dimension of input clusters."
    result = Tuple{Int, Int, Int}[]
    length(neighbors)==0 && return result
    for (i, indices) in enumerate(inrange(KDTree(convert(Matrix{Float}, cluster₂)), convert(Matrix{Float}, cluster₁), max(values(neighbors)...)+atol, true))
        for j in indices
            dist = zero(promote_type(eltype(cluster₁), eltype(cluster₂)))
            for k = 1:size(cluster₁, 1)
                dist = dist + (cluster₂[k, j]-cluster₁[k, i])^2
            end
            dist = sqrt(dist)
            for (nb, len) in neighbors
                if isapprox(len, dist, atol=atol, rtol=rtol)
                    push!(result, (nb, i, j))
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
    AbstractPID <: OperatorUnit

Abstract point id.
"""
abstract type AbstractPID <: OperatorUnit end

"""
    PID <: AbstractPID

Point id.
"""
struct PID <: AbstractPID
    site::Int
end

"""
    CPID{S} <: AbstractPID

Composite point id.
"""
struct CPID{S} <: AbstractPID
    scope::S
    site::Int
end

"""
    CPID(site::Int)
    CPID(; scope='T', site::Int=1)

Construct a composite point id.
"""
@inline CPID(site::Int) = CPID('T', site)
@inline CPID(; scope='T', site::Int=1) = CPID(scope, site)

"""
    AbstractBond{N, P<:AbstractPID, D<:Number, R}

Abstract bond.
"""
abstract type AbstractBond{N, P<:AbstractPID, D<:Number, R} end
@inline Base.:(==)(b₁::AbstractBond, b₂::AbstractBond) = ==(efficientoperations, b₁, b₂)
@inline Base.isequal(b₁::AbstractBond, b₂::AbstractBond) = isequal(efficientoperations, b₁, b₂)

"""
    length(bond::AbstractBond) -> Int
    length(::Type{<:AbstractBond{N, <:AbstractPID, <:Number, R} where N}) where R -> Int

Get the number of points of a bond.
"""
@inline Base.length(bond::AbstractBond) = length(typeof(bond))
@inline Base.length(::Type{<:AbstractBond{N, <:AbstractPID, <:Number, R} where N}) where R = R

"""
    eltype(bond::AbstractBond)
    eltype(::Type{<:AbstractBond{N, P}}) where {N, P<:AbstractPID}

Get the type of the points contained in a bond.
"""
@inline Base.eltype(bond::AbstractBond) = eltype(typeof(bond))
@inline Base.eltype(::Type{<:AbstractBond{N, P, D}}) where {N, P<:AbstractPID, D<:Number} = Point{N, P, D}

"""
    iterate(bond::AbstractBond, state=1)

Iterate over the points in a bond.
"""
@inline Base.iterate(bond::AbstractBond, state=1) = state<=length(bond) ? (bond[state], state+1) : nothing

"""
    dimension(bond::AbstractBond) -> Int
    dimension(::Type{<:AbstractBond{N}}) where N -> Int

Get the space dimension of a concrete bond.
"""
@inline dimension(bond::AbstractBond) = dimension(typeof(bond))
@inline dimension(::Type{<:AbstractBond{N}}) where N = N

"""
    pidtype(bond::AbstractBond)
    pidtype(::Type{<:AbstractBond{N, P} where N}) where {P<:AbstractPID}

Get the pid type of a concrete bond.
"""
@inline pidtype(bond::AbstractBond) = pidtype(typeof(bond))
@inline pidtype(::Type{<:AbstractBond{N, P} where N}) where {P<:AbstractPID} = P

"""
    dtype(bond::AbstractBond)
    dtype(::Type{<:AbstractBond{N, <:AbstractPID, D} where N}) where {D<:Number}

Get the data type of a bond.
"""
@inline dtype(bond::AbstractBond) = dtype(typeof(bond))
@inline dtype(::Type{<:AbstractBond{N, <:AbstractPID, D} where N}) where {D<:Number} = D

"""
    rank(bond::AbstractBond) -> Int
    rank(::Type{<:AbstractBond{N, <:AbstractPID, <:Number, R} where N}) where R -> Int

Get the rank of a bond.
"""
@inline rank(bond::AbstractBond) = rank(typeof(bond))
@inline rank(::Type{<:AbstractBond{N, <:AbstractPID, <:Number, R} where N}) where R = R

"""
    Point{N, P<:AbstractPID, D<:Number} <: AbstractBond{N, P, D, 1}

Labeled point.
"""
struct Point{N, P<:AbstractPID, D<:Number} <: AbstractBond{N, P, D, 1}
    pid::P
    rcoord::SVector{N, D}
    icoord::SVector{N, D}
    Point(pid::AbstractPID, rcoord::SVector{N, D}, icoord::SVector{N, D}) where {N, D<:Number} = new{N, typeof(pid), D}(pid, rcoord, icoord)
end
Base.show(io::IO, p::Point) = @printf io "Point(%s, %s, %s)" p.pid p.rcoord p.icoord

"""
    Point(pid::AbstractPID, rcoord::SVector{N, D}, icoord::SVector{N, D}) where {N, D<:Number}
    Point(pid::AbstractPID, rcoord::NTuple{N, <:Number}, icoord::NTuple{N, <:Number}=ntuple(i->0, N)) where N
    Point(pid::AbstractPID, rcoord::AbstractVector{<:Number}, icoord::AbstractVector{<:Number}=zero(SVector{length(rcoord), Int}))
    Point{N}(pid::AbstractPID, rcoord::AbstractVector{<:Number}, icoord::AbstractVector{<:Number}=zero(SVector{N, Int})) where N

Construct a labeled point.
"""
@inline function Point(pid::AbstractPID, rcoord::NTuple{N, <:Number}, icoord::NTuple{N, <:Number}=ntuple(i->0, N)) where N
    datatype = promote_type(eltype(rcoord), eltype(icoord))
    return Point(pid, convert(SVector{N, datatype}, rcoord), convert(SVector{N, datatype}, icoord))
end
@inline function Point(pid::AbstractPID, rcoord::AbstractVector{<:Number}, icoord::AbstractVector{<:Number}=zero(SVector{length(rcoord), Int}))
    return Point{length(rcoord)}(pid, rcoord, icoord)
end
@inline function Point{N}(pid::AbstractPID, rcoord::AbstractVector{<:Number}, icoord::AbstractVector{<:Number}=zero(SVector{N, Int})) where N
    @assert length(rcoord)==length(icoord)==N "Point error: mismatched length of input rcoord and icoord."
    datatype = promote_type(eltype(rcoord), eltype(icoord))
    return Point(pid, convert(SVector{N, datatype}, rcoord), convert(SVector{N, datatype}, icoord))
end

"""
    getindex(p::Point, i::Integer) -> Point

Get the ith point of a point, which by definition i should only be 1 and the result is the point itself.
"""
@inline function Base.getindex(p::Point, i::Integer)
    @assert i==1 "getindex error: for Point, the index should only be 1."
    return p
end

"""
    kind(::Point) -> 0
    kind(::Type{<:Point}) -> 0

Get the bond kind of a point, which is defined to be 0.
"""
@inline kind(::Point) = 0
@inline kind(::Type{<:Point}) = 0

"""
    isintracell(p::Point) -> Bool

Judge whether a point is intra the unitcell.
"""
@inline isintracell(p::Point) = isapprox(norm(p.icoord), 0.0, atol=atol, rtol=rtol)

"""
    Bond{N, P<:AbstractPID, D<:Number} <: AbstractBond{N, P, D, 2}

A bond in a lattice.
"""
struct Bond{N, P<:AbstractPID, D<:Number} <: AbstractBond{N, P, D, 2}
    neighbor::Int
    spoint::Point{N, P, D}
    epoint::Point{N, P, D}
end
Base.show(io::IO, bond::Bond) = @printf io "Bond(%s, %s, %s)" bond.neighbor bond.spoint bond.epoint

"""
    reverse(bond::Bond) -> Bond

Get the reversed bond.
"""
@inline Base.reverse(bond::Bond) = Bond(bond.neighbor, bond.epoint, bond.spoint)

"""
    getindex(bond::Bond, i::Integer) -> Point

Get the ith point of a bond.
"""
@inline function Base.getindex(bond::Bond, i::Integer)
    @assert i==1 || i==2 "getindex error: for Bond, the index should only be 1 or 2."
    i==1 && return bond.epoint
    return bond.spoint
end

"""
    kind(bond::Bond) -> Int

Get the bond kind of a bond.
"""
@inline kind(bond::Bond) = bond.neighbor

"""
    rcoord(bond::Bond) -> SVector

Get the rcoord of the bond.
"""
@inline rcoord(bond::Bond) = bond.epoint.rcoord - bond.spoint.rcoord

"""
    icoord(bond::Bond) -> SVector

Get the icoord of the bond.
"""
@inline icoord(bond::Bond) = bond.epoint.icoord - bond.spoint.icoord

"""
    isintracell(bond::Bond) -> Bool

Judge whether a bond is intra the unit cell of a lattice.
"""
@inline isintracell(bond::Bond) = isapprox(bond|>icoord|>norm, 0.0, atol=atol, rtol=rtol)

"""
    AbstractLattice{N, P<:AbstractPID, D<:Number}

Abstract type for all lattices.

It should have the following contents:
- `name::Symbol`: the name of the lattice
- `pids::Vector{P}`: the pids of the lattice
- `rcoords::Matrix{D}`: the rcoords of the lattice
- `icoords::Matrix{D}`: the icoords of the lattice
- `vectors::Vector{SVector{N, D}}`: the translation vectors of the lattice
- `reciprocals::Vector{SVector{N, D}}`: the reciprocals of the lattice
- `neighbors::Dict{Int, Float}`: the order-distance map of the nearest neighbors of the lattice
"""
abstract type AbstractLattice{N, P<:AbstractPID, D<:Number} end
@inline contentnames(::Type{<:AbstractLattice}) = (:name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)
@inline Base.:(==)(lattice1::AbstractLattice, lattice2::AbstractLattice) = ==(efficientoperations, lattice1, lattice2)
@inline Base.isequal(lattice1::AbstractLattice, lattice2::AbstractLattice) = isequal(efficientoperations, lattice1, lattice2)
function Base.show(io::IO, lattice::AbstractLattice)
    @printf io "%s(%s)\n" lattice|>typeof|>nameof getcontent(lattice, :name)
    len = length(lattice)
    if len > 0
        @printf io "  with %s %s:\n" len len==1 ? "point" : "points"
        for i = 1:len
            @printf io "    %s\n" lattice[LatticeIndex{'P'}(i)]
        end
    end
    len = length(getcontent(lattice, :vectors))
    if len > 0
        @printf io "  with %s translation %s:\n" len len==1 ? "vector" : "vectors"
        for i = 1:len
            @printf io "    %s\n" getcontent(lattice, :vectors)[i]
        end
    end
    len = length(getcontent(lattice, :neighbors))
    if len > 0
        @printf io "  with %s %s of nearest neighbors:\n" len len==1 ? "order" : "orders"
        for (order, distance) in getcontent(lattice, :neighbors)
            @printf io "    %s => %s\n" order distance
        end
    end
end

"""
    length(lattice::AbstractLattice) -> Int

Get the number of points contained in a lattice.
"""
@inline Base.length(lattice::AbstractLattice) = length(getcontent(lattice, :pids))

"""
    keytype(lattice::AbstractLattice)
    keytype(::Type{<:AbstractLattice{N, P} where N}) where {P<:AbstractPID}

Get the pid type of the lattice.
"""
@inline Base.keytype(lattice::AbstractLattice) = keytype(typeof(lattice))
@inline Base.keytype(::Type{<:AbstractLattice{N, P} where N}) where {P<:AbstractPID} = P

"""
    valtype(lattice::AbstractLattice)
    valtype(::Type{<:AbstractLattice{N, P}}) where {N, P<:AbstractPID}

Get the point type of the lattice.
"""
@inline Base.valtype(lattice::AbstractLattice) = valtype(typeof(lattice))
@inline Base.valtype(::Type{<:AbstractLattice{N, P, D}}) where {N, P<:AbstractPID, D<:Number} = Point{N, P, D}

"""
    dimension(lattice::AbstractLattice) -> Int
    dimension(::Type{<:AbstractLattice{N}}) where N -> Int

Get the space dimension of the lattice.
"""
@inline dimension(lattice::AbstractLattice) = dimension(typeof(lattice))
@inline dimension(::Type{<:AbstractLattice{N}}) where {N} = N

"""
    dtype(lattice::AbstractLattice)
    dtype(::Type{<:AbstractLattice{N, <:AbstractPID, D} where N}) where {D<:Number}

Get the data type of the coordinates of a lattice.
"""
@inline dtype(lattice::AbstractLattice) = dtype(typeof(lattice))
@inline dtype(::Type{<:AbstractLattice{N, <:AbstractPID, D} where N}) where {D<:Number} = D

"""
    nneighbor(lattice::AbstractLattice) -> Int

Get the highest order of nearest neighbors.
"""
@inline nneighbor(lattice::AbstractLattice) = max(keys(getcontent(lattice, :neighbors))...)

"""
    LatticeIndex{Kind}(index::Union{AbstractPID, Int}) where Kind

Lattice index.

`Kind` must be one of the followings:
* 'R': for getting the rcoord of a lattice
* 'I': for getting the icoord of a lattice
* 'P': for getting the point of a lattice
"""
struct LatticeIndex{Kind, I<:Union{AbstractPID, Int}}
    index::I
    function LatticeIndex{Kind}(index::Union{AbstractPID, Int}) where Kind
        @assert Kind in ('R', 'I', 'P') "LatticeIndex error: wrong input Kind($Kind)."
        new{Kind, typeof(index)}(index)
    end
end

"""
    getindex(lattice::AbstractLattice, i::LatticeIndex{'R'}) -> SVector
    getindex(lattice::AbstractLattice, i::LatticeIndex{'I'}) -> SVector
    getindex(lattice::AbstractLattice, i::LatticeIndex{'P'}) -> Point

Get a rcoord, an icoord or a point of a lattice according to the type of the input index.
"""
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'R', Int}) = latticestaticcoords(getcontent(lattice, :rcoords), i.index, lattice|>dimension|>Val)
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'I', Int}) = latticestaticcoords(getcontent(lattice, :icoords), i.index, lattice|>dimension|>Val)
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'P', Int}) = Point(getcontent(lattice, :pids)[i.index], lattice[LatticeIndex{'R'}(i.index)], lattice[LatticeIndex{'I'}(i.index)])
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'R', <:AbstractPID}) = lattice[LatticeIndex{'R'}(findfirst(isequal(i.index), getcontent(lattice, :pids)))]
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'I', <:AbstractPID}) = lattice[LatticeIndex{'I'}(findfirst(isequal(i.index), getcontent(lattice, :pids)))]
@inline Base.getindex(lattice::AbstractLattice, i::LatticeIndex{'P', <:AbstractPID}) = lattice[LatticeIndex{'P'}(findfirst(isequal(i.index), getcontent(lattice, :pids)))]
@inline latticestaticcoords(coords::Matrix{D}, i::Int, ::Val{N}) where {D<:Number, N} = SVector{N, D}(ntuple(j->coords[j, i], Val(N)))

abstract type LatticeBonds{R} end
@inline Base.eltype(::Type{<:LatticeBonds{0}}) = AbstractBond
@inline Base.eltype(::Type{<:LatticeBonds{1}}) = Point
@inline Base.eltype(::Type{<:LatticeBonds{2}}) = Bond
@inline Base.eltype(::Type{L}, ::B) where {L<:AbstractLattice, B<:LatticeBonds} = eltype(B){L|>dimension, L|>keytype, L|>dtype}
@inline Base.eltype(::Type{L}, ::Val{B}) where {L<:AbstractLattice, B} = eltype(typeof(B)){L|>dimension, L|>keytype, L|>dtype}
struct AllBonds <: LatticeBonds{0} end
struct ZerothBonds <: LatticeBonds{1} end
struct InsideBonds <: LatticeBonds{2} end
struct AcrossBonds <: LatticeBonds{2} end
"""
    allbonds

Indicate that all bonds are inquired.
"""
const allbonds = AllBonds()
"""
    zerothbonds

Indicate that zeroth bonds, i.e. the points are inquired.
"""
const zerothbonds = ZerothBonds()
"""
    insidebonds

Indicate that bonds inside the unitcell are inquired, which do not contain those across the periodic boundaries.
"""
const insidebonds = InsideBonds()
"""
    acrossbonds

Indicate that bonds across the unitcell are inquired, which are in fact those across the periodic boundaries.
"""
const acrossbonds = AcrossBonds()

"""
    latticebondsstructure(::Type{<:AbstractLattice}) -> SimpleTree{LatticeBonds, Nothing}

The tree structure of the lattice bonds.
"""
@generated function latticebondsstructure(::Type{<:AbstractLattice})
    structure = SimpleTree{LatticeBonds, Nothing}()
    push!(structure, allbonds, nothing)
    push!(structure, allbonds, zerothbonds, nothing)
    push!(structure, allbonds, insidebonds, nothing)
    push!(structure, allbonds, acrossbonds, nothing)
    return structure
end

"""
    expand(::Type{L}, ::Val{VS}) where {L<:AbstractLattice, VS} -> Tuple

Expand the lattice bond types to the leaf level.
"""
@generated function expand(::Type{L}, ::Val{VS}) where {L<:AbstractLattice, VS}
    leaves = LatticeBonds[]
    structure = latticebondsstructure(L)
    for V in VS
        for latticebonds in keys(structure, simpletreedepth, V)
            isleaf(structure, latticebonds) && push!(leaves, latticebonds)
        end
    end
    return Tuple(leaves)
end

"""
    bonds(lattice::AbstractLattice, inquiry=allbonds) -> Vector{eltype(lattice|>typeof, inquiry)}
    bonds(lattice::AbstractLattice, inquiries...) -> Vector{mapreduce(inquiry->eltype(lattice|>typeof, inquiry), typejoin, inquiries)}

Generate the required bonds of a lattice.
"""
@inline bonds(lattice::AbstractLattice) = bonds(lattice, allbonds)
@inline bonds(lattice::AbstractLattice, inquiry) = bonds!(eltype(lattice|>typeof, inquiry)[], lattice, inquiry)
@inline bonds(lattice::AbstractLattice, inquiries...) = bonds!(mapreduce(inquiry->eltype(lattice|>typeof, inquiry), typejoin, inquiries)[], lattice, inquiries...)

"""
    bonds!(bonds::Vector, lattice::AbstractLattice, inquiries::LatticeBonds...) -> Vector
    bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{zerothbonds}) -> Vector
    bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{insidebonds}) -> Vector
    bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{acrossbonds}) -> Vector

Generate the required bonds of a lattice and append them to the input bonds.
"""
@inline bonds!(bonds::Vector, lattice::AbstractLattice) = bonds!(bonds, lattice, allbonds)
function bonds!(bonds::Vector, lattice::AbstractLattice, inquiries::LatticeBonds...)
    leaves = expand(typeof(lattice), Val(inquiries))
    for i = 1:length(leaves)
        bonds!(bonds, lattice, leaves[i]|>Val)
    end
    return bonds
end
function bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{zerothbonds})
    for i = 1:length(lattice)
        push!(bonds, lattice[LatticeIndex{'P'}(i)])
    end
    return bonds
end
function bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{insidebonds})
    for (neighbor, sindex, eindex) in interlinks(getcontent(lattice, :rcoords), getcontent(lattice, :rcoords), getcontent(lattice, :neighbors))
        if sindex < eindex
            spoint = lattice[LatticeIndex{'P'}(sindex)]
            epoint = lattice[LatticeIndex{'P'}(eindex)]
            push!(bonds, Bond(neighbor, spoint, epoint))
        end
    end
    return bonds
end
function bonds!(bonds::Vector, lattice::AbstractLattice, ::Val{acrossbonds})
    nnb = nneighbor(lattice)
    translations = reshape(product((-nnb:nnb for i = 1:length(getcontent(lattice, :vectors)))...)|>collect, :)
    for translation in translations
        remove = ntuple(i->-translation[i], length(translation))
        filter!(≠(remove), translations)
    end
    if length(translations) > 0
        dim = lattice |> dimension |> Val
        superrcoords = tile(getcontent(lattice, :rcoords), getcontent(lattice, :vectors), translations)
        supericoords = tile(getcontent(lattice, :icoords), getcontent(lattice, :vectors), translations)
        for (neighbor, sindex, eindex) in interlinks(getcontent(lattice, :rcoords), superrcoords, getcontent(lattice, :neighbors))
            spoint = Point(getcontent(lattice, :pids)[sindex],
                        latticestaticcoords(getcontent(lattice, :rcoords), sindex, dim),
                        latticestaticcoords(getcontent(lattice, :icoords), sindex, dim)
                        )
            epoint = Point(getcontent(lattice, :pids)[(eindex-1)%length(lattice)+1],
                        latticestaticcoords(superrcoords, eindex, dim),
                        latticestaticcoords(supericoords, eindex, dim)
                        )
            push!(bonds, Bond(neighbor, spoint, epoint))
        end
    end
    return bonds
end

"""
    @recipe plot(lattice::AbstractLattice, bondtype::LatticeBonds=allbonds, filter::Function=bond->true)

Define the recipe for the visualization of a lattice.
"""
@recipe function plot(lattice::AbstractLattice, bondtype::LatticeBonds=allbonds, filter::Function=bond->true)
    title := String(getcontent(lattice, :name))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    for leafbond in expand(typeof(lattice), Val((bondtype,)))
        if leafbond == zerothbonds
            @series begin
                seriestype := :scatter
                rcoords = NTuple{dimension(lattice), dtype(lattice)}[]
                for i = 1:length(lattice)
                    point = lattice[LatticeIndex{'P'}(i)]
                    filter(point) && push!(rcoords, Tuple(point.rcoord))
                end
                rcoords
            end
        else
            @series begin
                data = Vector{Float64}[]
                for i = 1:dimension(lattice)
                    push!(data, Float64[])
                end
                for bond in bonds(lattice, leafbond)
                    if filter(bond)
                        for i = 1:dimension(lattice)
                            for point in bond
                                push!(data[i], point.rcoord[i])
                            end
                            push!(data[i], NaN)
                        end
                    end
                end
                Tuple(data)
            end
        end
    end
end

"""
    Lattice{N, P<:AbstractPID, D<:Number} <: AbstractLattice{N, P, D}

Simplest lattice.

A simplest lattice can be constructed from its contents, i.e. pids, rcoords and icoords, or from a couple of points, or from a couple of sublattices.
"""
struct Lattice{N, P<:AbstractPID, D<:Number} <: AbstractLattice{N, P, D}
    name::Symbol
    pids::Vector{P}
    rcoords::Matrix{D}
    icoords::Matrix{D}
    vectors::Vector{SVector{N, D}}
    reciprocals::Vector{SVector{N, D}}
    neighbors::Dict{Int, Float}
    function Lattice{N}(name::Symbol,
            pids::Vector{<:AbstractPID}, rcoords::AbstractMatrix{<:Number}, icoords::AbstractMatrix{<:Number},
            vectors::AbstractVector{<:AbstractVector{<:Number}},
            neighbors::Union{Dict{Int, <:Real}, Int}=1;
            coordination::Int=8
            ) where N
        @assert N==size(rcoords, 1)==size(icoords, 1) && length(pids)==size(rcoords, 2)==size(icoords, 2) "Lattice error: shape mismatched."
        isa(neighbors, Int) && (neighbors = Dict(i=>minlen for (i, minlen) in enumerate(minimumlengths(rcoords, vectors, neighbors, coordination=coordination))))
        datatype = promote_type(Float, eltype(rcoords), eltype(icoords), eltype(eltype(vectors)))
        rcoords = convert(Matrix{datatype}, rcoords)
        icoords = convert(Matrix{datatype}, icoords)
        vectors = convert(Vector{SVector{N, datatype}}, vectors)
        recipls = convert(Vector{SVector{N, datatype}}, reciprocals(vectors))
        new{N, eltype(pids), datatype}(name, pids, rcoords, icoords, vectors, recipls, neighbors)
    end
end

"""
    Lattice{N}(name::Symbol,
        pids::Vector{<:AbstractPID}, rcoords::AbstractMatrix{<:Number}, icoords::AbstractMatrix{<:Number},
        vectors::AbstractVector{<:AbstractVector{<:Number}},
        neighbors::Union{Dict{Int, <:Real}, Int}=1;
        coordination::Int=8
        ) where N
    Lattice(name::Symbol,
        points::AbstractVector{<:Point};
        vectors::AbstractVector{<:AbstractVector{<:Number}}=SVector{0, SVector{points|>eltype|>dimension, points|>eltype|>dtype}}(),
        neighbors::Union{Dict{Int, <:Real}, Int}=1,
        coordination::Int=8
        )
    Lattice(name::Symbol,
        sublattices::AbstractVector{<:AbstractLattice};
        vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0, SVector{sublattices|>eltype|>dimension, sublattices|>eltype|>dtype}}(),
        neighbors::Union{Dict{Int, <:Real}, Int}=1,
        coordination::Int=8
        )

Construct a lattice.
"""
function Lattice(name::Symbol,
        points::AbstractVector{<:Point};
        vectors::AbstractVector{<:AbstractVector{<:Number}}=SVector{0, SVector{points|>eltype|>dimension, points|>eltype|>dtype}}(),
        neighbors::Union{Dict{Int, <:Real}, Int}=1,
        coordination::Int=8
        )
    datatype = dtype(eltype(points))
    pids = Vector{points|>eltype|>pidtype}(undef, points|>length)
    rcoords = zeros(datatype, points|>eltype|>dimension, points|>length)
    icoords = zeros(datatype, points|>eltype|>dimension, points|>length)
    for i = 1:length(points)
        pids[i] = points[i].pid
        rcoords[:, i] = points[i].rcoord
        icoords[:, i] = points[i].icoord
    end
    return Lattice{points|>eltype|>dimension}(name, pids, rcoords, icoords, vectors, neighbors, coordination=coordination)
end
function Lattice(name::Symbol,
        sublattices::AbstractVector{<:AbstractLattice};
        vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0, SVector{sublattices|>eltype|>dimension, sublattices|>eltype|>dtype}}(),
        neighbors::Union{Dict{Int, <:Real}, Int}=1,
        coordination::Int=8
        )
    if length(sublattices) > 1
        @assert mapreduce(sl->dimension(sl)==dimension(sublattices[1]), &, sublattices) "Lattice error: all sublattices should have the same dimension."
        @assert mapreduce(sl->keytype(sl)===keytype(sublattices[1]), &, sublattices) "Lattice error: all sublattices should have the same keytype."
    end
    datatype = dtype(eltype(sublattices))
    len = mapreduce(length, +, sublattices)
    pids = Vector{sublattices|>eltype|>keytype}(undef, len)
    rcoords = zeros(datatype, sublattices|>eltype|>dimension, len)
    icoords = zeros(datatype, sublattices|>eltype|>dimension, len)
    pos = 1
    for sublattice in sublattices
        pids[pos:pos+length(sublattice)-1] = sublattice.pids
        rcoords[:, pos:pos+length(sublattice)-1] = sublattice.rcoords
        icoords[:, pos:pos+length(sublattice)-1] = sublattice.icoords
        pos += length(sublattice)
    end
    Lattice{sublattices|>eltype|>dimension}(name, pids, rcoords, icoords, vectors, neighbors, coordination=coordination)
end

"""
    Lattice(lattice::Lattice, translations::Translations; pidmap::Function=pidmap, icoordmap::Function=icoordmap, coordination::Int=8)

Construct a lattice from the translations of another.
"""
function Lattice(lattice::Lattice, translations::Translations; pidmap::Function=pidmap, icoordmap::Function=icoordmap, coordination::Int=8)
    name = Symbol(@sprintf "%s(%s)" lattice.name translations)
    rcoords = tile(lattice.rcoords, lattice.vectors, translations)
    vectors = SVector{dimension(lattice), dtype(lattice)}[]
    for (i, vector) in enumerate(lattice.vectors)
        translations.boundaries[i]=='P' && push!(vectors, vector*translations.ranges[i])
    end
    len = length(lattice)*length(translations)
    pids = Vector{keytype(lattice)}(undef, len)
    icoords = zeros(dtype(lattice), dimension(lattice), len)
    count = 1
    for translation in translations
        for index = 1:length(lattice)
            pids[count] = pidmap(lattice, index, translations, translation)
            icoords[:, count] .= icoordmap(lattice, index, translations, translation)
            count += 1
        end
    end
    return Lattice{dimension(lattice)}(name, pids, rcoords, icoords, vectors, lattice.neighbors, coordination=coordination)
end
function pidmap(lattice::Lattice, index::Int, translations::Translations{N}, translation::NTuple{N, Int}) where N
    pid = lattice.pids[index]
    sequence = findfirst(translation, translations)
    return replace(pid, site=pid.site+(sequence-1)*length(lattice))
end
icoordmap(lattice::Lattice, index::Int, ::Translations{N}, ::NTuple{N, Int}) where N = lattice[LatticeIndex{'I'}(index)]

"""
    SuperLattice(name::Symbol,
        sublattices::AbstractVector{<:AbstractLattice};
        vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0, SVector{sublattices|>eltype|>dimension, sublattices|>eltype|>dtype}}(),
        neighbors::Dict{Int, <:Real}=Dict{Int, Float}()
        )

SuperLattice that is composed of several sublattices.
"""
struct SuperLattice{L<:AbstractLattice, N, P<:AbstractPID, D<:Number} <: AbstractLattice{N, P, D}
    sublattices::Vector{L}
    name::Symbol
    pids::Vector{P}
    rcoords::Matrix{D}
    icoords::Matrix{D}
    vectors::Vector{SVector{N, D}}
    reciprocals::Vector{SVector{N, D}}
    neighbors::Dict{Int, Float}
    function SuperLattice(name::Symbol,
            sublattices::AbstractVector{<:AbstractLattice};
            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0, SVector{sublattices|>eltype|>dimension, sublattices|>eltype|>dtype}}(),
            neighbors::Dict{Int, <:Real}=Dict{Int, Float}()
            )
        if length(sublattices) > 1
            @assert mapreduce(sl->dimension(sl)==dimension(sublattices[1]), &, sublattices) "SuperLattice error: all sublattices should have the same dimension."
            @assert mapreduce(sl->keytype(sl)===keytype(sublattices[1]), &, sublattices) "SuperLattice error: all sublattices should have the same keytype."
        end
        @assert all(length(sublattice.vectors) == 0 for sublattice in sublattices) "SuperLattice error: all sublattices should assume open boundaries."
        datatype = sublattices|>eltype|>dtype
        len = mapreduce(length, +, sublattices)
        pids = Vector{sublattices|>eltype|>keytype}(undef, len)
        rcoords = zeros(datatype, sublattices|>eltype|>dimension, len)
        icoords = zeros(datatype, sublattices|>eltype|>dimension, len)
        pos = 1
        for sublattice in sublattices
            pids[pos:pos+length(sublattice)-1] = sublattice.pids
            rcoords[:, pos:pos+length(sublattice)-1] = sublattice.rcoords
            icoords[:, pos:pos+length(sublattice)-1] = sublattice.icoords
            pos += length(sublattice)
        end
        vectors = convert(Vector{SVector{sublattices|>eltype|>dimension, datatype}}, vectors)
        recipls = convert(Vector{SVector{sublattices|>eltype|>dimension, datatype}}, reciprocals(vectors))
        new{sublattices|>eltype, sublattices|>eltype|>dimension, sublattices|>eltype|>keytype, datatype}(sublattices, name, pids, rcoords, icoords, vectors, recipls, neighbors)
    end
end
@inline contentnames(::Type{<:SuperLattice}) = (:sublattices, :name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)

"""
    latticetype(sl::SuperLattice)
    latticetype(::Type{<:SuperLattice{L}}) where {L<:AbstractLattice}

Get the sublattice type of a superlattice.
"""
@inline latticetype(sl::SuperLattice) = latticetype(typeof(sl))
@inline latticetype(::Type{<:SuperLattice{L}}) where {L<:AbstractLattice} = L

struct IntraBonds <: LatticeBonds{2} end
struct InterBonds <: LatticeBonds{2} end
"""
    intrabonds

Indicate that bonds intra the sublattices are inquired.
!!! note
    These bonds do not contain those across the periodic boundaries.
"""
const intrabonds = IntraBonds()
"""
    interbonds

Indicate that bonds inter the sublattices are inquired.
!!! note
    These bonds do not contain those across the periodic boundaries.
"""
const interbonds = InterBonds()

"""
    latticebondsstructure(::Type{<:SuperLattice}) -> SimpleTree{LatticeBonds, Nothing}

The tree structure of the lattice bonds.
"""
@generated function latticebondsstructure(::Type{<:SuperLattice})
    structure = SimpleTree{LatticeBonds, Nothing}()
    push!(structure, allbonds, nothing)
    push!(structure, allbonds, zerothbonds, nothing)
    push!(structure, allbonds, insidebonds, nothing)
    push!(structure, allbonds, acrossbonds, nothing)
    push!(structure, insidebonds, intrabonds, nothing)
    push!(structure, insidebonds, interbonds, nothing)
    return structure
end

"""
    bonds!(bonds::Vector, lattice::SuperLattice, ::Val{intrabonds}) -> Vector
    bonds!(bonds::Vector, lattice::SuperLattice, ::Val{interbonds}) -> Vector

Generate the required bonds of a superlattice and append them to the input bonds.
"""
function bonds!(bonds::Vector, lattice::SuperLattice, ::Val{intrabonds})
    for sublattice in lattice.sublattices
        bonds!(bonds, sublattice, insidebonds)
    end
    return bonds
end
function bonds!(bonds::Vector, lattice::SuperLattice, ::Val{interbonds})
    dim = lattice |> dimension |> Val
    for (i, j) in Combinations{2}(1:length(lattice.sublattices))
        sub1, sub2 = lattice.sublattices[i], lattice.sublattices[j]
        for (neighbor, sindex, eindex) in interlinks(sub1.rcoords, sub2.rcoords, lattice.neighbors)
            spoint = Point(sub1.pids[sindex], latticestaticcoords(sub1.rcoords, sindex, dim), latticestaticcoords(sub1.icoords, sindex, dim))
            epoint = Point(sub2.pids[eindex], latticestaticcoords(sub2.rcoords, eindex, dim), latticestaticcoords(sub2.icoords, eindex, dim))
            push!(bonds, Bond(neighbor, spoint, epoint))
        end
    end
    return bonds
end

"""
    Cylinder{P}(name::Symbol, block::AbstractMatrix{<:Real}, translation::SVector{N, <:Real};
        vector::Union{AbstractVector{<:Real}, Nothing}=nothing,
        neighbors::Union{Dict{Int, <:Real}, Int}=1
        ) where {P<:CPID, N}

Cylinder of 1d and quasi 2d lattices.
"""
mutable struct Cylinder{P<:CPID, N, D<:Number} <: AbstractLattice{N, P, D}
    block::Matrix{D}
    translation::SVector{N, D}
    name::Symbol
    pids::Vector{P}
    rcoords::Matrix{D}
    icoords::Matrix{D}
    vectors::Vector{SVector{N, D}}
    reciprocals::Vector{SVector{N, D}}
    neighbors::Dict{Int, Float}
    function Cylinder{P}(name::Symbol, block::AbstractMatrix{<:Number}, translation::SVector{N, <:Number};
            vector::AbstractVector{<:Number}=SVector{N, Float}(),
            neighbors::Union{Dict{Int, <:Real}, Int}=1
            ) where {P<:CPID, N}
        pids = Vector{P}[]
        datatype = promote_type(Float, eltype(block), eltype(translation), eltype(vector))
        rcoords = zeros(datatype, N, 0)
        icoords = zeros(datatype, N, 0)
        vectors = [convert(SVector{N, datatype}, vector)]
        recipls = convert(Vector{SVector{N, datatype}}, reciprocals(vectors))
        isa(neighbors, Int) && (neighbors = Dict{Int, Float}(i=>Inf for i = 1:neighbors))
        new{P, N, datatype}(block, translation, name, pids, rcoords, icoords, vectors, recipls, neighbors)
    end
end
@inline contentnames(::Type{<:Cylinder}) = (:block, :translation, :name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)

"""
    insert!(cylinder::Cylinder, ps::S...; cut=length(cylinder)÷2+1, scopes::Union{<:AbstractVector{S}, Nothing}=nothing, coordination=8) where S -> Cylinder

Insert a couple of blocks into a cylinder.

The position of the cut of the cylinder is specified by the keyword argument `cut`, which is the center of the cylinder by default. All pids corresponding to a same newly inserted block share the same scope, which is specified by the parameter `ps`. Optionally, the scopes of the old pids in the cylinder can be replaced if the parameter `scopes` is assigned other than `nothing`. Note the length of `ps` is equal to the number of newly inserted blocks, while that of `scopes` should be equal to the old length of the cylinder.
"""
function Base.insert!(cylinder::Cylinder, ps::S...; cut=length(cylinder)÷2+1, scopes::Union{<:AbstractVector{S}, Nothing}=nothing, coordination=8) where S
    @assert S<:fieldtype(cylinder|>keytype, :scope) "insert! error: mismatched type of input scopes and old scopes."
    @assert 1<=cut<=length(cylinder)+1 "insert! error: wrong cut($cut), which should be in [1, $(length(cylinder)+1)]."
    @assert length(cylinder)%size(cylinder.block, 2)==0 "insert! error: wrong cylinder length."
    isnothing(scopes) || @assert length(scopes)==length(cylinder) "insert! error: mismatched length of input scopes and cylinder."
    blcklen = size(cylinder.block, 2)
    blcknum = length(cylinder)÷blcklen + length(ps)
    rsltlen = blcklen * blcknum
    pids = keytype(cylinder)[]
    rcoords = zeros(dtype(cylinder), cylinder|>dimension, rsltlen)
    icoords = zeros(dtype(cylinder), cylinder|>dimension, rsltlen)
    for i = 1:cut-1
        push!(pids, isnothing(scopes) ? cylinder.pids[i] : replace(cylinder.pids[i], scope=scopes[i]))
    end
    for (i, p) in enumerate(ps)
        for j = 1:blcklen
            push!(pids, CPID(p, cut==length(cylinder)+1 ? j : (cylinder.pids[cut].site+j-2)%blcklen+1))
        end
    end
    for i = cut:length(cylinder)
        push!(pids, isnothing(scopes) ? cylinder.pids[i] : replace(cylinder.pids[i], scope=scopes[i]))
    end
    for i = 1:blcknum
        trannum = i - (blcknum+1)/2
        for j = 1:blcklen
            for k = 1:dimension(cylinder)
                rcoords[k, (i-1)*blcklen+j] = cylinder.block[k, j] + trannum*cylinder.translation[k]
            end
        end
    end
    cylinder.pids = pids
    cylinder.rcoords = rcoords
    cylinder.icoords = icoords
    if length(cylinder.neighbors)>0 && mapreduce(isequal(Inf), |, values(cylinder.neighbors))
        minlens = minimumlengths(rcoords, cylinder.vectors, cylinder|>nneighbor, coordination=coordination)
        cylinder.neighbors = Dict{Int, Float}(i=>len for (i, len) in enumerate(minlens))
    end
    return cylinder
end

"""
    (cylinder::Cylinder)(scopes::Any...; coordination::Int=8) -> Lattice

Construct a lattice from a cylinder with the assigned scopes.
"""
function (cylinder::Cylinder)(scopes::Any...; coordination::Int=8)
    @assert fieldtype(cylinder|>keytype, :scope)==scopes|>eltype "cylinder call error: wrong scope type."
    name = Symbol(@sprintf "%s%s" cylinder.name length(scopes))
    pids = keytype(cylinder)[CPID(scope, i) for scope in scopes for i = 1:size(cylinder.block, 2)]
    rcoords = tile(cylinder.block, [cylinder.translation], [(i,) for i = -(length(scopes)-1)/2:(length(scopes)-1)/2])
    icoords = zero(rcoords)
    vectors = cylinder.vectors
    neighbors = cylinder.neighbors
    if length(neighbors)>0 && mapreduce(isequal(Inf), |, values(neighbors))
        minlens = minimumlengths(rcoords, vectors, cylinder|>nneighbor, coordination=coordination)
        neighbors = Dict{Int, Float}(i=>len for (i, len) in enumerate(minlens))
    end
    Lattice{cylinder|>dimension}(name, pids, rcoords, icoords, vectors, neighbors)
end

"""
    Bonds{T, L<:AbstractLattice, BS<:Tuple{Vararg{Vector{<:AbstractBond}}}, B<:AbstractBond} <: AbstractVector{B}

A set of lattice bonds.

`Bonds` itself is an `AbstractVector` of `AbstractBond`. The need for such a struct is to ensure the type stability during the iteration over a set of different concrete bonds. Although the default `iterate` function does not achieve this goal, users can get it with the generated function trick. Besides, it provides a high level of management of different categories of bonds based on the `LatticeBonds` system.
"""
struct Bonds{T, L<:AbstractLattice, BS<:Tuple{Vararg{Vector{<:AbstractBond}}}, B<:AbstractBond} <: AbstractVector{B}
    bonds::BS
    function Bonds{T, L}(bonds::Tuple{Vararg{Vector{<:AbstractBond}}}) where {T, L<:AbstractLattice}
        @assert isa(T, Tuple{Vararg{LatticeBonds}}) && length(T)==length(bonds) "Bonds error: mismatched input types and bonds."
        new{T, L, typeof(bonds), mapreduce(eltype, typejoin, bonds)}(bonds)
    end
end
@inline Base.size(bs::Bonds) = (mapreduce(length, +, bs.bonds),)
@inline Base.:(==)(bonds1::Bonds, bonds2::Bonds) = ==(efficientoperations, bonds1, bonds2)
@inline Base.isequal(bonds1::Bonds, bonds2::Bonds) = isequal(efficientoperations, bonds1, bonds2)
@inline Base.summary(io::IO, bs::Bonds) = @printf io "%s-element Bonds" length(bs)
function Base.iterate(bs::Bonds, state=(1, 0))
    s1, s2 = state[1], state[2]
    while s1 <= length(bs.bonds)
        s1, s2 = (s2+1 > length(bs.bonds[s1])) ? (s1+1, 1) : (s1, s2+1)
        (s1>length(bs.bonds) || s2<=length(bs.bonds[s1])) && break
    end
    s1>length(bs.bonds) && return nothing
    return bs.bonds[s1][s2], (s1, s2)
end
function Base.getindex(bs::Bonds, i::Int)
    r, k = 1, i
    while r<=rank(bs) && k>length(bs.bonds[r])
        k = k - length(bs.bonds[r])
        r = r + 1
    end
    r>rank(bs) && error("getindex error: attempt to access $(length(bs))-element Bonds at index [$k].")
    return bs.bonds[r][k]
end

"""
    Bonds{T, L}(bonds::Tuple{Vararg{Vector{<:AbstractBond}}}) where {T, L<:AbstractLattice}
    Bonds(lattice::AbstractLattice, types::LatticeBonds...)

Construct a set of lattice bonds.
"""
@inline Bonds(lattice::AbstractLattice) = Bonds(lattice, allbonds)
@inline Bonds(lattice::AbstractLattice, types::LatticeBonds...) = Bonds(lattice, Val(types))
@generated function Bonds(lattice::AbstractLattice, ::Val{VS}) where VS
    types, bonds = [], []
    for V in expand(lattice, Val(VS))
        push!(types, V)
        push!(bonds, :(bonds(lattice, Val($V))))
    end
    types = Tuple(types)
    bonds = Expr(:tuple, bonds...)
    return :(Bonds{$types, typeof(lattice)}($bonds))
end

"""
    bondtypes(bs::Bonds) -> Tuple{Vararg{LatticeBonds}}
    bondtypes(::Type{<:Bonds{T}}) where T -> Tuple{Vararg{LatticeBonds}}

Get the bondtypes of a set of lattice bonds.
"""
@inline bondtypes(bs::Bonds) = bondtypes(typeof(bs))
@inline bondtypes(::Type{<:Bonds{T}}) where {T} = T

"""
    latticetype(bs::Bonds)
    latticetype(::Type{<:Bonds{T, L} where T}) where {L<:AbstractLattice}
"""
@inline latticetype(bs::Bonds) = latticetype(typeof(bs))
@inline latticetype(::Type{<:Bonds{T, L} where T}) where {L<:AbstractLattice} = L

"""
    rank(bs::Bonds) -> Int
    rank(::Type{<:Bonds{T}}) where T -> Int

Get the rank of a set of lattice bonds.
"""
@inline rank(bs::Bonds) = rank(typeof(bs))
@inline rank(::Type{<:Bonds{T}}) where T = length(T)

"""
    filter(lbs::LatticeBonds, bs::Bonds, choice::Union{Val{:include}, Val{:exclude}}=Val(:include)) -> Bonds
    filter(lbs::Tuple{Vararg{LatticeBonds}}, bs::Bonds, choice::Union{Val{:include}, Val{:exclude}}=Val(:include)) -> Bonds

Get a subset of a set of lattice bonds.

When `choice=Val(:include)`, the lattice bonds indicated by `lbs` will be selected; 
When `choice=Val(:exclude)`, the lattice bonds not indicated by `lbs` will be selected.
"""
@inline Base.filter(lbs::LatticeBonds, bs::Bonds, choice::Union{Val{:include}, Val{:exclude}}=Val(:include)) = bondfilter(Val((lbs,)), bs, choice)
@inline Base.filter(lbs::Tuple{Vararg{LatticeBonds}}, bs::Bonds, choice::Union{Val{:include}, Val{:exclude}} = Val(:include)) = bondfilter(Val(lbs), bs, choice)
@generated function bondfilter(::Val{VS}, bs::Bonds, ::Val{:include}) where VS
    types, bonds = [], []
    for V in expand(latticetype(bs), Val(VS))
        index = findfirst(isequal(V), bondtypes(bs))
        if isa(index, Int)
            push!(types, V)
            push!(bonds, :(bs.bonds[$index]))
        end
    end
    types = Tuple(types)
    bonds = Expr(:tuple, bonds...)
    return :(Bonds{$types, latticetype(bs)}($bonds))
end
@generated function bondfilter(::Val{VS}, bs::Bonds, ::Val{:exclude}) where VS
    types, bonds = [], []
    excludes = expand(latticetype(bs), Val(VS))
    for (i, V) in enumerate(bondtypes(bs))
        if V ∉ excludes
            push!(types, V)
            push!(bonds, :(bs.bonds[$i]))
        end
    end
    types = Tuple(types)
    bonds = Expr(:tuple, bonds...)
    return :(Bonds{$types, latticetype(bs)}($bonds))
end

"""
    filter(select::Function, bs::Bonds) -> Bonds

Get a filtered set of bonds by a select function.
"""
@inline Base.filter(select::Function, bs::Bonds) = Bonds{bondtypes(bs), latticetype(bs)}(map(bonds->filter(select, bonds), bs.bonds))

"""
    empty!(bs::Bonds) -> Bonds

Empty a set of lattice bonds.
"""
@inline Base.empty!(bs::Bonds) = (map(empty!, bs.bonds); bs)

"""
    empty(bs::Bonds) -> Bonds

Get an empty copy of a set of lattice bonds.
"""
@inline Base.empty(bs::Bonds) = Bonds{bondtypes(bs), latticetype(bs)}(map(empty, bs.bonds))

"""
    reset!(bs::Bonds, lattice::AbstractLattice) -> Bonds

Reset a set of lattice bonds by a new lattice.
"""
function reset!(bs::Bonds, lattice::AbstractLattice)
    @assert latticetype(bs) >: typeof(lattice) "reset! error: mismatched bonds and lattice."
    empty!(bs)
    map((bonds, bondtype)->bonds!(bonds, lattice, bondtype), bs.bonds, bondtypes(bs))
    return bs
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
        @assert isa(K, Symbol) || isa(K, Tuple{Symbol}) "BrillouinZone error: K must be a Symbol."
        reciprocals = vectorconvert(reciprocals)
        new{K, eltype(momenta), eltype(reciprocals)}(reciprocals, momenta)
    end
end
@inline contentnames(::Type{<:BrillouinZone}) = (:reciprocals, :content)
@inline getcontent(bz::BrillouinZone, ::Val{:content}) = (bz.momenta,)

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
    volume::V
    function ReciprocalZone{K}(momenta::AbstractVector, volume::Number) where K
        @assert isa(K, Symbol) || isa(K, Tuple{Symbol}) "ReciprocalZone error: K must be a Symbol."
        momenta = vectorconvert(momenta)
        new{K, eltype(momenta), typeof(volume)}(momenta, volume)
    end
end
@inline contentnames(::Type{<:ReciprocalZone}) = (:content, :volume)
@inline getcontent(rz::ReciprocalZone, ::Val{:content}) = (rz.momenta,)

"""
    ReciprocalZone(momenta::AbstractVector, volume::Number)
    ReciprocalZone(reciprocals::AbstractVector; length::Int=100)
    ReciprocalZone(reciprocals::AbstractVector, bounds::Segment...)
    ReciprocalZone(reciprocals::AbstractVector, bounds::Tuple{Vararg{Segment}})

    ReciprocalZone{K}(momenta::AbstractVector, volume::Number) where K
    ReciprocalZone{K}(reciprocals::AbstractVector; length::Int=100) where K
    ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Segment...) where K
    ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Tuple{Vararg{Segment}}) where K

Construct a rectangular zone in the reciprocal space.
"""
@inline ReciprocalZone(momenta::AbstractVector, volume::Number) = ReciprocalZone{:k}(momenta, volume)
@inline ReciprocalZone(reciprocals::AbstractVector; length::Int=100) = ReciprocalZone{:k}(vectorconvert(reciprocals); length=length)
@inline ReciprocalZone(reciprocals::AbstractVector, bounds::Segment...) = ReciprocalZone{:k}(vectorconvert(reciprocals), bounds)
@inline ReciprocalZone(reciprocals::AbstractVector, bounds::Tuple{Vararg{Segment}}) = ReciprocalZone{:k}(vectorconvert(reciprocals), bounds)
@inline function ReciprocalZone{K}(reciprocals::AbstractVector; length::Int=100) where K
    return ReciprocalZone{K}(vectorconvert(reciprocals), ntuple(i->Segment(-1//2, 1//2, length, ends=(true, false)), Base.length(reciprocals)))
end
@inline ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Segment...) where K = ReciprocalZone{K}(vectorconvert(reciprocals), bounds)
function ReciprocalZone{K}(reciprocals::AbstractVector, bounds::Tuple{Vararg{Segment}}) where K
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
    return ReciprocalZone{K}(momenta, ratio*volume(reciprocals))
end

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
        @assert isa(K, Symbol) || isa(K, Tuple{Symbol}) "ReciprocalPath error: K must be a Symbol."
        momenta = vectorconvert(momenta)
        new{K, eltype(momenta)}(momenta)
    end
end
@inline contentnames(::Type{<:ReciprocalPath}) = (:content,)
@inline getcontent(rp::ReciprocalPath, ::Val{:content}) = (rp.momenta,)

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
