module Spatials

using LinearAlgebra: norm,dot,cross,inv
using StaticArrays: SVector,SMatrix
using Printf: @printf
using NearestNeighbors: KDTree,knn,inrange
using Base.Iterators: flatten,product
using ...Prerequisites: atol,rtol,Float
using ...Prerequisites.TypeTraits: efficientoperations
using ...Prerequisites.SimpleTrees: SimpleTree,simpletreedepth,isleaf
using ...Mathematics.Combinatorics: Combinations
using ...Mathematics.AlgebraOverFields: SimpleID

import ...Interfaces: decompose,rank,dimension,kind,expand,reset!

export distance,azimuthd,azimuth,polard,polar,volume
export isparallel,isonline,isintratriangle,issubordinate
export reciprocals,translate,rotate,tile,minimumlengths
export intralinks,interlinks
export PID,AbstractBond,pidtype,Point,Bond,rcoord,icoord,isintracell
export AbstractLattice,nneighbor,LatticeIndex,bonds!,bonds,LatticeBonds,Bonds,latticetype,bondtypes
export latticebondsstructure,allbonds,zerothbonds,insidebonds,acrossbonds,intrabonds,interbonds
export Lattice,SuperLattice,Cylinder

"""
    distance(p1::AbstractVector{<:Real},p2::AbstractVector{<:Real}) -> Float

Get the distance between two points.

!!! note
    Compared to `norm(p1-p2)`, this function avoids the memory allocation for `p1-p2`, thus is more efficient.
"""
function distance(p1::AbstractVector{<:Real},p2::AbstractVector{<:Real})
    @assert length(p1)==length(p2) "distance error: dismatched length of input vectors."
    result=0
    for i=1:length(p1)
        result=result+(p1[i]-p2[i])^2
    end
    sqrt(result)
end

"""
    azimuthd(v::AbstractVector{<:Real}) -> Float

Get the azimuth angle in degrees of a vector.
"""
function azimuthd(v::AbstractVector{<:Real})
    @assert length(v) ∈ (1,2,3) "azimuthd error: wrong dimensioned input vector."
    result=acosd(v[1]/(length(v)==3 ? sqrt(v[1]^2+v[2]^2) : norm(v)))
    length(v)>1 && v[2]<0 && (result=360-result)
    return result
end

"""
    azimuth(v::AbstractVector{<:Real}) -> Float

Get the azimuth angle in radians of a vector.
"""
function azimuth(v::AbstractVector{<:Real})
    @assert length(v) ∈ (1,2,3) "azimuth error: wrong dimensioned input vector."
    result=acos(v[1]/(length(v)==3 ? sqrt(v[1]^2+v[2]^2) : norm(v)))
    length(v)>1 && v[2]<0 && (result=2pi-result)
    return result
end

"""
    polard(v::AbstractVector{<:Real}) -> Float

Get the polar angle in degrees of a vector.
"""
function polard(v::AbstractVector{<:Real})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acosd(v[3]/norm(v))
end

"""
    polar(v::AbstractVector{<:Real}) -> Float

Get the polar angle in radians of a vector.
"""
function polar(v::AbstractVector{<:Real})
    @assert length(v)==3 "polard error: wrong dimensioned input vector."
    return acos(v[3]/norm(v))
end

"""
    volume(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real}) -> Real

Get the volume spanned by three vectors.
"""
function volume(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real})
    if length(v1) ∈ (1,2) || length(v2) ∈ (1,2) || length(v3) ∈ (1,2)
        result=zero(v1|>eltype)
    elseif length(v1)==3 && length(v2)==3 && length(v3)==3
        result=v1[1]*(v2[2]*v3[3]-v2[3]*v3[2])+v1[2]*(v2[3]*v3[1]-v2[1]*v3[3])+v1[3]*(v2[1]*v3[2]-v2[2]*v3[1])
    else
        error("volume error: wrong dimensioned input vectors.")
    end
    return result
end

"""
    isparallel(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real};atol::Real=atol,rtol::Real=rtol) -> Int

Judge whether two vectors are parallel to each other with the given tolerance, `0` for not parallel, `1` for parallel and `-1` for antiparallel.
"""
function isparallel(v1::AbstractVector{<:Real},v2::AbstractVector{<:Real};atol::Real=atol,rtol::Real=rtol)
    norm1,norm2=norm(v1),norm(v2)
    if isapprox(norm1,0.0,atol=atol,rtol=rtol) || isapprox(norm2,0.0,atol=atol,rtol=rtol)
        result=1
    elseif length(v1)==length(v2)
        temp=dot(v1,v2)/norm1/norm2
        result=isapprox(temp,1,atol=atol,rtol=rtol) ? 1 : isapprox(temp,-1,atol=atol,rtol=rtol) ? -1 : 0
    else
        error("isparallel error: shape dismatch of the input vectors.")
    end
    return result
end

"""
    isonline(p::AbstractVector{<:Real},p1::AbstractVector{<:Real},p2::AbstractVector{<:Real};ends::Tuple{Bool,Bool}=(true,true),atol::Real=atol,rtol::Real=rtol) -> Bool

Judge whether a point is on a line segment whose end points are `p1` and `p2` with the given tolerance. `ends` defines whether the line segment should contain its ends.
"""
function isonline(p::AbstractVector{<:Real},p1::AbstractVector{<:Real},p2::AbstractVector{<:Real};ends::Tuple{Bool,Bool}=(true,true),atol::Real=atol,rtol::Real=rtol)
    @assert length(p)==length(p1)==length(p2) "isonline error: shape dismatch of input point and line segment."
    d1,d2,d=distance(p,p1),distance(p,p2),distance(p1,p2)
    isapprox(d1,0.0,atol=atol,rtol=rtol) && return ends[1]
    isapprox(d2,0.0,atol=atol,rtol=rtol) && return ends[2]
    return isapprox(d1+d2,d,atol=atol,rtol=rtol)
end

"""
    decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real}) -> Tuple{Float}
    decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real},v2::AbstractVector{<:Real}) -> Tuple{Float,Float}
    decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real}) -> Tuple{Float,Float,Float}

Decompose a vector with respect to input basis vectors.
"""
function decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real})
    @assert length(v0)==length(v1) "decompose error: dismatched length of input vectors."
    n0,n1=norm(v0),norm(v1)
    sign=dot(v0,v1)/n0/n1
    @assert isapprox(abs(sign),1.0,atol=atol,rtol=rtol) "decompose error: insufficient basis vectors."
    return (sign*n0/n1,)
end
function decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real},v2::AbstractVector{<:Real})
    @assert length(v0)==length(v1)==length(v2) "decompose error: dismatched length of input vectors."
    @assert length(v0)==2 || length(v0)==3 "decompose error: unsupported dimension($(length(v0))) of input vectors."
    if length(v0)==2
        det=v1[1]*v2[2]-v1[2]*v2[1]
        x1=(v0[1]*v2[2]-v0[2]*v2[1])/det
        x2=(v0[2]*v1[1]-v0[1]*v1[2])/det
    else
        v3=SVector{3,Float}(v1[2]*v2[3]-v1[3]*v2[2],v1[3]*v2[1]-v1[1]*v2[3],v1[1]*v2[2]-v1[2]*v2[1])
        x1,x2,x3=decompose(v0,v1,v2,v3)
        @assert isapprox(x3,0.0,atol=atol,rtol=rtol) "decompose error: insufficient basis vectors."
    end
    return x1,x2
end
function decompose(v0::AbstractVector{<:Real},v1::AbstractVector{<:Real},v2::AbstractVector{<:Real},v3::AbstractVector{<:Real})
    @assert length(v0)==length(v1)==length(v2)==length(v3) "decompose error: dismatched length of input vectors."
    @assert length(v0)==3 "decompose error: unsupported dimension($(length(v0))) of input vectors."
    V=volume(v1,v2,v3)
    r1=(v2[2]*v3[3]/V-v2[3]*v3[2]/V,v2[3]*v3[1]/V-v2[1]*v3[3]/V,v2[1]*v3[2]/V-v2[2]*v3[1]/V)
    r2=(v3[2]*v1[3]/V-v3[3]*v1[2]/V,v3[3]*v1[1]/V-v3[1]*v1[3]/V,v3[1]*v1[2]/V-v3[2]*v1[1]/V)
    r3=(v1[2]*v2[3]/V-v1[3]*v2[2]/V,v1[3]*v2[1]/V-v1[1]*v2[3]/V,v1[1]*v2[2]/V-v1[2]*v2[1]/V)
    return dot(r1,v0),dot(r2,v0),dot(r3,v0)
end

"""
    isintratriangle(p::AbstractVector{<:Real},
                    p1::AbstractVector{<:Real},
                    p2::AbstractVector{<:Real},
                    p3::AbstractVector{<:Real};
                    vertexes::NTuple{3,Bool}=(true,true,true),
                    edges::NTuple{3,Bool}=(true,true,true),
                    atol::Real=atol,
                    rtol::Real=rtol
                    ) -> Bool

Judge whether a point belongs to the interior of a triangle whose vertexes are `p1`, 'p2' and `p3` with the give tolerance. `vertexes` and `edges` define whether the interior should contain the vertexes or edges, respectively.
!!! note
    1. The vertexes are in the order (p1,p2,p3) and the edges are in the order (p1p2,p2p3,p3p1).
    2. The edges do not contain the vertexes.
"""
function isintratriangle(p::AbstractVector{<:Real},
                p1::AbstractVector{<:Real},p2::AbstractVector{<:Real},p3::AbstractVector{<:Real};
                vertexes::NTuple{3,Bool}=(true,true,true),
                edges::NTuple{3,Bool}=(true,true,true),
                atol::Real=atol,
                rtol::Real=rtol
                )
    @assert length(p)==length(p1)==length(p2)==length(p3) "isintratriangle error: shape dismatch of input point and triangle."
    @assert length(p)==2 || length(p)==3 "isintratriangle error: unsupported dimension($(length(p))) of input points."
    x=length(p)==2 ? decompose(SVector(p[1]-p1[1],p[2]-p1[2]),SVector(p2[1]-p1[1],p2[2]-p1[2]),SVector(p3[1]-p1[1],p3[2]-p1[2])) :
                     decompose(SVector(p[1]-p1[1],p[2]-p1[2],p[3]-p1[3]),SVector(p2[1]-p1[1],p2[2]-p1[2],p2[3]-p1[3]),SVector(p3[1]-p1[1],p3[2]-p1[2],p3[3]-p1[3]))
    x1_approx_0,x2_approx_0=isapprox(x[1],0.0,atol=atol,rtol=rtol),isapprox(x[2],0.0,atol=atol,rtol=rtol)
    x1_approx_1,x2_approx_1=isapprox(x[1],1.0,atol=atol,rtol=rtol),isapprox(x[2],1.0,atol=atol,rtol=rtol)
    x12_approx_1=isapprox(x[1]+x[2],1.0,atol=atol,rtol=rtol)
    x1_approx_0 && x2_approx_0 && return vertexes[1]
    x1_approx_1 && x2_approx_0 && return vertexes[2]
    x1_approx_0 && x2_approx_1 && return vertexes[3]
    x2_approx_0 && 0<x[1]<1 && return edges[1]
    x12_approx_1 && 0<x[1]<1 && 0<x[2]<1 && return edges[2]
    x1_approx_0 && 0<x[2]<1 && return edges[3]
    return 0<x[1]<1 && 0<x[2]<1 && x[1]+x[2]<1
end

"""
    issubordinate(rcoord::AbstractVector{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}};atol::Real=atol,rtol::Real=rtol) -> Bool

Judge whether a coordinate belongs to a lattice defined by `vectors` with the given tolerance.
"""
function issubordinate(rcoord::AbstractVector{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}};atol::Real=atol,rtol::Real=rtol)
    @assert length(rcoord) ∈ (1,2,3) "issubordinate error: only 1, 2 and 3 dimensional coordinates are supported."
    @assert length(vectors) ∈ (1,2,3) "issubordinate error: the number of input basis vectors must be 1, 2 or 3."
    fapprox=xi->isapprox(round(xi),xi,atol=atol,rtol=rtol)
    if length(vectors)==1
        result=mapreduce(fapprox,&,decompose(rcoord,vectors[1]))
    elseif length(vectors)==2
        result=mapreduce(fapprox,&,decompose(rcoord,vectors[1],vectors[2]))
    else
        result=mapreduce(fapprox,&,decompose(rcoord,vectors[1],vectors[2],vectors[3]))
    end
    return result
end

"""
    reciprocals(vectors::AbstractVector{AbstractVector{<:Real}}) -> Vector{Vector{Float}}

Get the reciprocals dual to the input vectors.
"""
function reciprocals(vectors::AbstractVector{<:AbstractVector{<:Real}})
    @assert length(vectors)<4 "reciprocals error: the number of input vectors should not be greater than 3."
    @assert mapreduce(v->length(v) in (1,2,3),&,vectors,init=true) "reciprocals error: all input vectors must be 1, 2 or 3 dimensional."
    result=Vector{Float}[]
    if length(vectors)==1
        push!(result,2pi/mapreduce(vi->vi^2,+,vectors[1])*vectors[1])
    elseif length(vectors)==2
        v1,v2=vectors[1],vectors[2]
        @assert length(v1)==length(v2) "reciprocals error: dismatched length of input vectors."
        if length(v1)==2
            det=2pi/(v1[1]*v2[2]-v1[2]*v2[1])
            push!(result,[det*v2[2],-det*v2[1]])
            push!(result,[-det*v1[2],det*v1[1]])
        else
            v3=SVector{3,Float}(v1[2]*v2[3]-v1[3]*v2[2],v1[3]*v2[1]-v1[1]*v2[3],v1[1]*v2[2]-v1[2]*v2[1])
            V=2pi/volume(v1,v2,v3)
            push!(result,[v2[2]*v3[3]*V-v2[3]*v3[2]*V,v2[3]*v3[1]*V-v2[1]*v3[3]*V,v2[1]*v3[2]*V-v2[2]*v3[1]*V])
            push!(result,[v3[2]*v1[3]*V-v3[3]*v1[2]*V,v3[3]*v1[1]*V-v3[1]*v1[3]*V,v3[1]*v1[2]*V-v3[2]*v1[1]*V])
        end
    elseif length(vectors)==3
        v1,v2,v3=vectors
        V=2pi/volume(v1,v2,v3)
        push!(result,[v2[2]*v3[3]*V-v2[3]*v3[2]*V,v2[3]*v3[1]*V-v2[1]*v3[3]*V,v2[1]*v3[2]*V-v2[2]*v3[1]*V])
        push!(result,[v3[2]*v1[3]*V-v3[3]*v1[2]*V,v3[3]*v1[1]*V-v3[1]*v1[3]*V,v3[1]*v1[2]*V-v3[2]*v1[1]*V])
        push!(result,[v1[2]*v2[3]*V-v1[3]*v2[2]*V,v1[3]*v2[1]*V-v1[1]*v2[3]*V,v1[1]*v2[2]*V-v1[2]*v2[1]*V])
    end
    return result
end

"""
    translate(cluster::AbstractMatrix{<:Real},vector::AbstractVector{<:Real}) -> Matrix{vector|>eltype}

Get the translated cluster of the original one by a vector.
"""
translate(cluster::AbstractMatrix{<:Real},vector::AbstractVector{<:Real})=cluster.+reshape(vector,(vector|>length,1))

"""
    rotate(cluster::AbstractMatrix{<:Real},angle::Real;axis::Tuple{Union{AbstractVector{<:Real},Nothing},Tuple{<:Real,<:Real}}=(nothing,(0,0))) -> Matrix{Float}

Get a rotated cluster of the original one by a certain angle around an axis.

The axis is determined by a point it gets through (`nothing` can be used to denote the origin), and its polar as well as azimuth angles in radians. The default axis is the z axis.
!!! note
    1. The result is given by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
    2. Only 2 and 3 dimensional vectors can be rotated.
    3. When the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.
"""
function rotate(cluster::AbstractMatrix{<:Real},angle::Real;axis::Tuple{Union{AbstractVector{<:Real},Nothing},Tuple{<:Real,<:Real}}=(nothing,(0,0)))
    @assert size(cluster,1) in (2,3) "rotate error: only 2 and 3 dimensional vectors can be rotated."
    center,theta,phi=(axis[1]===nothing ? zeros(size(cluster,1)) : axis[1]),axis[2][1],axis[2][2]
    @assert length(center)==size(cluster,1) "rotate error: dismatched shape of the input cluster and the point on axis."
    length(center)==2 && @assert isapprox(theta,0,atol=atol) && isapprox(phi,0,atol=atol) "rotate error: both the polar and azimuth of the axis for 2d vectors must be 0."
    cosθ,sinθ=cos(angle),sin(angle)
    k,w=SVector{3,Float}(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)),zeros(3)
    result=zeros(size(cluster))
    for i=1:size(cluster,2)
        for j=1:size(cluster,1)
            w[j]=cluster[j,i]-center[j]
        end
        inner,outer=dot(k,w),(k[2]*w[3]-k[3]*w[2],k[3]*w[1]-k[1]*w[3],k[1]*w[2]-k[2]*w[1])
        for j=1:size(cluster,1)
            result[j,i]=w[j]*cosθ+outer[j]*sinθ+k[j]*inner*(1-cosθ)+center[j]
        end
    end
    return result
end

"""
    tile(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},translations::NTuple{M,NTuple{N,<:Real}}=()) where {N,M} -> Matrix{Float}

Tile a supercluster by translations of the input cluster.

Basically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by `vectors` and each set of the translation indices contained in `translations`. When translation vectors are empty, a copy of the original cluster will be returned.
"""
function tile(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},translations::NTuple{M,NTuple{N,<:Real}}=()) where {N,M}
    N==0 && return copy(cluster)
    @assert length(vectors)==N "tile error: dismatched shape of input vectors and translations."
    supercluster=zeros(size(cluster,1),size(cluster,2)*length(translations))
    disp=zeros(size(cluster,1))
    for (i,translation) in enumerate(translations)
        for i=1:length(disp)
            disp[i]=0.0
            for j=1:length(vectors)
                disp[i]+=vectors[j][i]*translation[j]
            end
        end
        for j=1:size(cluster,2)
            col=(i-1)*size(cluster,2)+j
            for row=1:size(cluster,1)
                supercluster[row,col]=cluster[row,j]+disp[row]
            end
        end
    end
    return supercluster
end

"""
    minimumlengths(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},nneighbor::Int=1;coordination::Int=8) -> Vector{Float}

Use kdtree to search the lowest several minimum bond lengths within a lattice translated by a cluster.

When the translation vectors are not empty, the lattice will be considered periodic in the corresponding directions. Otherwise the lattice will be open in all directions. To search for the bonds accorss the periodic boundaries, the cluster will be pretranslated to become a supercluster, which has open boundaries but is large enough to contain all the nearest neighbors within the required order. The `coordination` parameter sets the average number of each order of nearest neighbors. If it is to small, larger bond lengths may not be searched, and the result will contain `Inf`. This is a sign that you may need a larger `coordination`. Another situation that `Inf` appears in the result occurs when the minimum lengths are searched in open lattices. Indeed, the cluster may be too small so that the required order just goes beyond it. In this case the warning message can be safely ignored.
"""
function minimumlengths(cluster::AbstractMatrix{<:Real},vectors::AbstractVector{<:AbstractVector{<:Real}},nneighbor::Int=1;coordination::Int=8)
    @assert nneighbor>=0 "minimumlengths error: input nneighbor must be non negative."
    result=[Inf for i=1:nneighbor]
    if size(cluster,2)>0
        translations=reshape(product((-nneighbor:nneighbor for i=1:length(vectors))...)|>collect,:)
        for translation in translations
            if length(translation)>0 && mapreduce(x->x!=0,|,translation)
                remove=ntuple(i->-translation[i],length(translation))
                filter!(key->key!=remove,translations)
            end
        end
        supercluster=tile(cluster,vectors,Tuple(translations))
        for len in flatten(knn(KDTree(supercluster),cluster,nneighbor>0 ? min(nneighbor*coordination,size(supercluster,2)) : 1,true)[2])
            for (i,minlen) in enumerate(result)
                if isapprox(len,minlen,atol=atol)
                    break
                elseif 0.0<len<minlen
                    nneighbor>0 && (result[i+1:nneighbor]=result[i:nneighbor-1])
                    result[i]=len
                    break
                end
            end
        end
        nneighbor>0 && mapreduce(isequal(Inf),|,result) && @warn "minimumlengths warning: larger(>$coordination) coordination or smaller(<$nneighbor) nneighbor may be needed."
    end
    return result
end

"""
    intralinks( cluster::AbstractMatrix{<:Real},
                vectors::AbstractVector{<:AbstractVector{<:Real}},
                neighbors::Dict{Int,Float},
                maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))
                ) where N -> Vector{Tuple{Int,Int,Int,SubArray{<:Real,1}}}

Use kdtree to get the intracluster nearest neighbors.

As is similar to [`minimumlengths`](@ref), when `vectors` is nonempty, the cluster assumes periodic boundaries. `neighbors` provides the map between the bond length and the order of nearest neighbors. Note only those with the lengths present in `neighbors` will be included in the result. `maxtranslations` determines the maximum number of translations along those directions specified by `vectors` when the tiled supercluster is construted (See [`minimumlengths`](@ref) for the explanation of the method for periodic lattices).
"""
function intralinks(    cluster::AbstractMatrix{<:Real},
                        vectors::AbstractVector{<:AbstractVector{<:Real}},
                        neighbors::Dict{Int,Float},
                        maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))
                        ) where N
    @assert length(vectors)==N "intralinks error: dismatched shape of input vectors and maxtranslations."
    result=Tuple{Int,Int,Int,Vector{Float}}[]
    length(neighbors)==0 && return result
    translations=reshape(product((-nnb:nnb for nnb in maxtranslations)...)|>collect,:)
    for translation in translations
        if length(translation)>0 && mapreduce(x->x!=0,|,translation)
            remove=ntuple(i->-translation[i],length(translation))
            filter!(key->key!=remove,translations)
        end
    end
    sort!(translations,by=norm)
    translations=Tuple(translations)
    supercluster=tile(cluster,vectors,translations)
    disps=tile(zero(cluster),vectors,translations)
    for (i,indices) in enumerate(inrange(KDTree(supercluster),cluster,max(values(neighbors)...)+atol,true))
        for j in indices
            if i<j
                dist=0.0
                for k=1:size(cluster,1)
                    dist=dist+(supercluster[k,j]-cluster[k,i])^2
                end
                dist=sqrt(dist)
                for (nb,len) in neighbors
                    isapprox(len,dist,atol=atol) && (push!(result,(nb,i,(j-1)%size(cluster,2)+1,disps[:,j]));break)
                end
            end
        end
    end
    return result
end

"""
    interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float}) -> Vector{Tuple{Int,Int,Int}}

Use kdtree to get the intercluster nearest neighbors.
"""
function interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float})
    @assert size(cluster1,1)==size(cluster2,1) "interlinks error: dismatched space dimension of input clusters."
    result=Tuple{Int,Int,Int}[]
    length(neighbors)==0 && return result
    for (i,indices) in enumerate(inrange(KDTree(cluster2),cluster1,max(values(neighbors)...)+atol,true))
        for j in indices
            dist=0.0
            for k=1:size(cluster1,1)
                dist=dist+(cluster2[k,j]-cluster1[k,i])^2
            end
            dist=sqrt(dist)
            for (nb,len) in neighbors
                isapprox(len,dist,atol=atol) && (push!(result,(nb,i,j));break)
            end
        end
    end
    return result
end

"""
    PID(scope,site::Int)
    PID(site::Int)
    PID(;scope="tz",site::Int=1)

The id of a point.
"""
struct PID{S} <: SimpleID
    scope::S
    site::Int
end
PID(site::Int)=PID('T',site)
PID(;scope='T',site::Int=1)=PID(scope,site)
Base.fieldnames(pid::PID)=(:scope,:site)

"""
    AbstractBond{N,P<:PID,R}

Abstract bond.
"""
abstract type AbstractBond{N,P<:PID,R} end

"""
    ==(b1::AbstractBond{N,P,R},b2::AbstractBond{N,P,R}) where {N,P<:PID,R} -> Bool
    isequal(b1::AbstractBond{N,P,R},b2::AbstractBond{N,P,R}) where {N,P<:PID,R} -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(b1::AbstractBond{N,P,R},b2::AbstractBond{N,P,R}) where {N,P<:PID,R} = ==(efficientoperations,b1,b2)
Base.isequal(b1::AbstractBond{N,P,R},b2::AbstractBond{N,P,R}) where {N,P<:PID,R}=isequal(efficientoperations,b1,b2)

"""
    length(bond::AbstractBond) -> Int
    length(::Type{<:AbstractBond{N,<:PID,R} where N}) where R -> Int

Get the number of points of a bond.
"""
Base.length(bond::AbstractBond)=bond|>typeof|>length
Base.length(::Type{<:AbstractBond{N,<:PID,R} where N}) where R=R

"""
    eltype(bond::AbstractBond)
    eltype(::Type{<:AbstractBond{N,P}}) where {N,P<:PID}

Get the eltype of a bond.
"""
Base.eltype(bond::AbstractBond)=bond|>typeof|>eltype
Base.eltype(::Type{<:AbstractBond{N,P}}) where {N,P<:PID}=Point{N,P}

"""
    dimension(bond::AbstractBond) -> Int
    dimension(::Type{<:AbstractBond{N}}) where N -> Int

Get the space dimension of a concrete bond.
"""
dimension(bond::AbstractBond)=bond|>typeof|>dimension
dimension(::Type{<:AbstractBond{N}}) where N=N

"""
    pidtype(bond::AbstractBond)
    pidtype(::Type{<:AbstractBond{N,P} where N}) where {P<:PID}

Get the pid type of a concrete bond.
"""
pidtype(bond::AbstractBond)=bond|>typeof|>pidtype
pidtype(::Type{<:AbstractBond{N,P} where N}) where {P<:PID}=P

"""
    rank(bond::AbstractBond) -> Int
    rank(::Type{<:AbstractBond{N,<:PID,R} where N}) where R -> Int

Get the rank of a bond.
"""
rank(bond::AbstractBond)=bond|>typeof|>rank
rank(::Type{<:AbstractBond{N,<:PID,R} where N}) where R=R

"""
    Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N
    Point(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N
    Point(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))
    Point{N}(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{N,Float})) where N

Labeled point.
"""
struct Point{N,P<:PID} <: AbstractBond{N,P,1}
    pid::P
    rcoord::SVector{N,Float}
    icoord::SVector{N,Float}
    Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N=new{N,pid|>typeof}(pid,rcoord,icoord)
end
Point(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N=Point(pid,convert(SVector{N,Float},rcoord),convert(SVector{N,Float},icoord))
Point(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))=Point{length(rcoord)}(pid,rcoord,icoord)
function Point{N}(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{N,Float})) where N
    @assert length(rcoord)==length(icoord)==N "Point error: dismatched length of input rcoord and icoord."
    Point(pid,convert(SVector{N,Float},rcoord),convert(SVector{N,Float},icoord))
end

"""
    show(io::IO,p::Point)

Show a labeled point.
"""
Base.show(io::IO,p::Point)=@printf io "Point(%s,[%s],[%s])" p.pid join(string.(p.rcoord),",") join(string.(p.icoord),",")

"""
    iterate(p::Point,state=1)

Iterate over the point.
"""
Base.iterate(p::Point,state=1)=state==1 ? (p,state+1) : nothing

"""
    kind(::Point) -> 0
    kind(::Type{<:Point}) -> 0

Get the bond kind of a point, which is defined to be 0.
"""
kind(::Point)=0
kind(::Type{<:Point})=0

"""
    Bond(neighbor::Int,spoint::Point,epoint::Point)

A bond in a lattice.
"""
struct Bond{N,P<:PID} <: AbstractBond{N,P,2}
    neighbor::Int
    spoint::Point{N,P}
    epoint::Point{N,P}
end

"""
    show(io::IO,bond::Bond)

Show a bond.
"""
Base.show(io::IO,bond::Bond)=@printf io "Bond(%s,%s,%s)" bond.neighbor bond.spoint bond.epoint

"""
    reverse(bond::Bond) -> Bond

Get the reversed bond.
"""
Base.reverse(bond::Bond)=Bond(bond.neighbor,bond.epoint,bond.spoint)

"""
    iterate(bond::Bond,state=1)

Iterate over the points in a bond.
"""
Base.iterate(bond::Bond,state=1)=state==1 ? (bond.epoint,state+1) : state==2 ? (bond.spoint,state+1) : nothing

"""
    kind(bond::Bond) -> Int

Get the bond kind of a bond.
"""
kind(bond::Bond)=bond.neighbor

"""
    rcoord(bond::Bond) -> SVector

Get the rcoord of the bond.
"""
rcoord(bond::Bond)=bond.epoint.rcoord-bond.spoint.rcoord

"""
    icoord(bond::Bond) -> SVector

Get the icoord of the bond.
"""
icoord(bond::Bond)=bond.epoint.icoord-bond.spoint.icoord

"""
    isintracell(bond::Bond) -> Bool

Judge whether a bond is intra the unit cell of a lattice.
"""
isintracell(bond::Bond)=isapprox(bond|>icoord|>norm,0.0,atol=atol)

"""
    AbstractLattice{P<:PID,N}

Abstract type for all lattices.

It should have the following attributes
- `name::String`: the name of the lattice
- `pids::Vector{P}`: the pids of the lattice
- `rcoords::Matrix{Float}`: the rcoords of the lattice
- `icoords::Matrix{Float}`: the icoords of the lattice
- `vectors::Vector{SVector{N,Float}}`: the translation vectors of the lattice
- `reciprocals::Vector{SVector{N,Float}}`: the reciprocals of the lattice
- `neighbors::Dict{Int,Float}`: the order-distance map of the nearest neighbors of the lattice
"""
abstract type AbstractLattice{N,P<:PID} end

"""
    show(io::IO,lattice::AbstractLattice)

Show a lattice.
"""
function Base.show(io::IO,lattice::AbstractLattice)
    @printf io "%s(%s)\n" lattice|>typeof|>nameof lattice.name
    if length(lattice)>0
        @printf io "  with %s %s:\n" length(lattice) length(lattice)==1 ? "point" : "points"
        for i=1:length(lattice)
            @printf io "    %s\n" lattice[LatticeIndex{'P'}(i)]
        end
    end
    if length(lattice.vectors)>0
        @printf io "  with %s translation %s:\n" length(lattice.vectors) length(lattice.vectors)==1 ? "vector" : "vectors"
        for i=1:length(lattice.vectors)
            @printf io "    [%s]\n" join(lattice.vectors[i],",")
        end
    end
    if length(lattice.neighbors)>0
        @printf io "  with %s %s of nearest neighbors:\n" length(lattice.neighbors) length(lattice.neighbors)==1 ? "order" : "orders"
        for (order,distance) in lattice.neighbors
            @printf io "    %s=>%s\n" order distance
        end
    end
end

"""
    length(lattice::AbstractLattice) -> Int

Get the number of points contained in a lattice.
"""
Base.length(lattice::AbstractLattice)=length(lattice.pids)

"""
    ==(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool
    isequal(lattice1::AbstractLattice,lattice2::AbstractLattice) -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(lattice1::AbstractLattice,lattice2::AbstractLattice) = ==(efficientoperations,lattice1,lattice2)
Base.isequal(lattice1::AbstractLattice,lattice2::AbstractLattice)=isequal(efficientoperations,lattice1,lattice2)

"""
    keytype(lattice::AbstractLattice)
    keytype(::Type{<:AbstractLattice{N,P} where N}) where {P<:PID}

Get the pid type of the lattice.
"""
Base.keytype(lattice::AbstractLattice)=lattice|>typeof|>keytype
Base.keytype(::Type{<:AbstractLattice{N,P} where N}) where {P<:PID}=P

"""
    valtype(lattice::AbstractLattice)
    valtype(::Type{<:AbstractLattice{N,P}}) where {N,P<:PID}

Get the point type of the lattice.
"""
Base.valtype(lattice::AbstractLattice)=lattice|>typeof|>valtype
Base.valtype(::Type{<:AbstractLattice{N,P}}) where {N,P<:PID}=Point{N,P}

"""
    dimension(lattice::AbstractLattice) -> Int
    dimension(::Type{<:AbstractLattice{N}}) where N -> Int

Get the space dimension of the lattice.
"""
dimension(lattice::AbstractLattice)=lattice|>typeof|>dimension
dimension(::Type{<:AbstractLattice{N}}) where N=N

"""
    nneighbor(lattice::AbstractLattice) -> Int

Get the highest order of nearest neighbors.
"""
nneighbor(lattice::AbstractLattice)=max(keys(lattice.neighbors)...)

"""
    LatticeIndex{Kind}(index::Union{PID,Int}) where Kind

Lattice index.

`Kind` must be one of the followings:
* 'R': for getting the rcoord of a lattice
* 'I': for getting the icoord of a lattice
* 'P': for getting the point of a lattice
"""
struct LatticeIndex{Kind,I<:Union{PID,Int}}
    index::I
    function LatticeIndex{Kind}(index::Union{PID,Int}) where Kind
        @assert Kind in ('R','I','P') "LatticeIndex error: wrong input Kind($Kind)."
        new{Kind,typeof(index)}(index)
    end
end

"""
    getindex(lattice::AbstractLattice,i::LatticeIndex{'R'}) -> SVector
    getindex(lattice::AbstractLattice,i::LatticeIndex{'I'}) -> SVector
    getindex(lattice::AbstractLattice,i::LatticeIndex{'P'}) -> Point

Get a rcoord, an icoord or a point of a lattice according to the type of the input index.
"""
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'R',Int})=latticestaticcoords(lattice.rcoords,i.index,lattice|>dimension|>Val)
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'I',Int})=latticestaticcoords(lattice.icoords,i.index,lattice|>dimension|>Val)
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'P',Int})=Point(lattice.pids[i.index],lattice[LatticeIndex{'R'}(i.index)],lattice[LatticeIndex{'I'}(i.index)])
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'R',<:PID})=lattice[LatticeIndex{'R'}(findfirst(isequal(i.index),lattice.pids))]
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'I',<:PID})=lattice[LatticeIndex{'I'}(findfirst(isequal(i.index),lattice.pids))]
Base.getindex(lattice::AbstractLattice,i::LatticeIndex{'P',<:PID})=lattice[LatticeIndex{'P'}(findfirst(isequal(i.index),lattice.pids))]
@generated function latticestaticcoords(coords::Matrix{Float},i::Int,::Val{N}) where N
    exprs=[:(coords[$j,i]) for j=1:N]
    return :(SVector{N,Float}($(exprs...)))
end

abstract type LatticeBonds{R} end
Base.eltype(::Type{<:LatticeBonds{0}})=AbstractBond
Base.eltype(::Type{<:LatticeBonds{1}})=Point
Base.eltype(::Type{<:LatticeBonds{2}})=Bond
Base.eltype(::Type{L},::B) where {L<:AbstractLattice,B<:LatticeBonds}=eltype(B){L|>dimension,L|>keytype}
Base.eltype(::Type{L},::Val{B}) where {L<:AbstractLattice,B}=eltype(typeof(B)){L|>dimension,L|>keytype}
struct AllBonds <: LatticeBonds{0} end
struct ZerothBonds <: LatticeBonds{1} end
struct InsideBonds <: LatticeBonds{2} end
struct AcrossBonds <: LatticeBonds{2} end
"""
    allbonds

Indicate that all bonds are inquired.
"""
const allbonds=AllBonds()
"""
    zerothbonds

Indicate that zeroth bonds, i.e. the points are inquired.
"""
const zerothbonds=ZerothBonds()
"""
    insidebonds

Indicate that bonds inside the unitcell are inquired, which do not contain those across the periodic boundaries.
"""
const insidebonds=InsideBonds()
"""
    acrossbonds

Indicate that bonds across the unitcell are inquired, which are in fact those across the periodic boundaries.
"""
const acrossbonds=AcrossBonds()

"""
    latticebondsstructure(::Type{<:AbstractLattice}) -> SimpleTree{LatticeBonds,Nothing}

The tree structure of the lattice bonds.
"""
@generated function latticebondsstructure(::Type{<:AbstractLattice})
    structure=SimpleTree{LatticeBonds,Nothing}()
    push!(structure,allbonds,nothing)
    push!(structure,allbonds,zerothbonds,nothing)
    push!(structure,allbonds,insidebonds,nothing)
    push!(structure,allbonds,acrossbonds,nothing)
    return structure
end

"""
    expand(::Type{L},::Val{VS}) where {L<:AbstractLattice,VS} -> Tuple

Expand the lattice bond types to the leaf level.
"""
@generated function expand(::Type{L},::Val{VS}) where {L<:AbstractLattice,VS}
    leaves=LatticeBonds[]
    structure=latticebondsstructure(L)
    for V in VS
        for latticebonds in keys(structure,simpletreedepth,V)
            isleaf(structure,latticebonds) && push!(leaves,latticebonds)
        end
    end
    return Tuple(leaves)
end

"""
    bonds(lattice::AbstractLattice,inquiry=allbonds) -> Vector{eltype(lattice|>typeof,inquiry)}
    bonds(lattice::AbstractLattice,inquiries...) -> Vector{mapreduce(inquiry->eltype(lattice|>typeof,inquiry),typejoin,inquiries)}

Generate the required bonds of a lattice.
"""
bonds(lattice::AbstractLattice)=bonds(lattice,allbonds)
bonds(lattice::AbstractLattice,inquiry)=bonds!(eltype(lattice|>typeof,inquiry)[],lattice,inquiry)
bonds(lattice::AbstractLattice,inquiries...)=bonds!(mapreduce(inquiry->eltype(lattice|>typeof,inquiry),typejoin,inquiries)[],lattice,inquiries...)

"""
    bonds!(bonds::Vector,lattice::AbstractLattice,inquiries::LatticeBonds...) -> Vector
    bonds!(bonds::Vector,lattice::AbstractLattice,::Val{zerothbonds}) -> Vector
    bonds!(bonds::Vector,lattice::AbstractLattice,::Val{insidebonds}) -> Vector
    bonds!(bonds::Vector,lattice::AbstractLattice,::Val{acrossbonds}) -> Vector

Generate the required bonds of a lattice and append them to the input bonds.
"""
bonds!(bonds::Vector,lattice::AbstractLattice)=bonds!(bonds,lattice,allbonds)
function bonds!(bonds::Vector,lattice::AbstractLattice,inquiries::LatticeBonds...)
    leaves=expand(typeof(lattice),Val(inquiries))
    for i=1:length(leaves)
        bonds!(bonds,lattice,leaves[i]|>Val)
    end
    return bonds
end
function bonds!(bonds::Vector,lattice::AbstractLattice,::Val{zerothbonds})
    for i=1:length(lattice)
        push!(bonds,lattice[LatticeIndex{'P'}(i)])
    end
    return bonds
end
function bonds!(bonds::Vector,lattice::AbstractLattice,::Val{insidebonds})
    for (neighbor,sindex,eindex) in interlinks(lattice.rcoords,lattice.rcoords,lattice.neighbors)
        if sindex<eindex
            spoint=lattice[LatticeIndex{'P'}(sindex)]
            epoint=lattice[LatticeIndex{'P'}(eindex)]
            push!(bonds,Bond(neighbor,spoint,epoint))
        end
    end
    return bonds
end
function bonds!(bonds::Vector,lattice::AbstractLattice,::Val{acrossbonds})
    nnb=lattice|>nneighbor
    translations=reshape(product((-nnb:nnb for i=1:length(lattice.vectors))...)|>collect,:)
    for translation in translations
        remove=ntuple(i->-translation[i],length(translation))
        filter!(key->key!=remove,translations)
    end
    if length(translations)>0
        dim=lattice|>dimension|>Val
        superrcoords=tile(lattice.rcoords,lattice.vectors,Tuple(translations))
        supericoords=tile(lattice.icoords,lattice.vectors,Tuple(translations))
        for (neighbor,sindex,eindex) in interlinks(lattice.rcoords,superrcoords,lattice.neighbors)
            spoint=Point(lattice.pids[sindex],latticestaticcoords(lattice.rcoords,sindex,dim),latticestaticcoords(lattice.icoords,sindex,dim))
            epoint=Point(lattice.pids[(eindex-1)%size(lattice.rcoords,2)+1],latticestaticcoords(superrcoords,eindex,dim),latticestaticcoords(supericoords,eindex,dim))
            push!(bonds,Bond(neighbor,spoint,epoint))
        end
    end
    return bonds
end

"""
    Lattice{N}( name::String,
                pids::Vector{<:PID},
                rcoords::AbstractMatrix{<:Real},
                icoords::AbstractMatrix{<:Real},
                vectors::AbstractVector{<:AbstractVector{<:Real}},
                neighbors::Union{Dict{Int,<:Real},Int}=1;
                coordination::Int=8
                ) where N
    Lattice(    name::String,
                points::AbstractVector{<:Point};
                vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{points|>eltype|>dimension,Float}}(),
                neighbors::Union{Dict{Int,<:Real},Int}=1,
                coordination::Int=8
                )
    Lattice(    name::String,
                sublattices::AbstractVector{<:AbstractLattice};
                vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                neighbors::Union{Dict{Int,<:Real},Int}=1,
                coordination::Int=8
                )

Simplest lattice.

A simplest lattice can be construted from its contents, i.e. pids, rcoords and icoords, or from a couple of points, or from a couple of sublattices.
"""
struct Lattice{N,P<:PID} <: AbstractLattice{N,P}
    name::String
    pids::Vector{P}
    rcoords::Matrix{Float}
    icoords::Matrix{Float}
    vectors::Vector{SVector{N,Float}}
    reciprocals::Vector{SVector{N,Float}}
    neighbors::Dict{Int,Float}
    function Lattice{N}(name::String,
                        pids::Vector{<:PID},
                        rcoords::AbstractMatrix{<:Real},
                        icoords::AbstractMatrix{<:Real},
                        vectors::AbstractVector{<:AbstractVector{<:Real}},
                        neighbors::Union{Dict{Int,<:Real},Int}=1;
                        coordination::Int=8
                ) where N
        @assert N==size(rcoords,1)==size(icoords,1) && length(pids)==size(rcoords,2)==size(icoords,2) "Lattice error: dismatched shape of input pids, rcoords and icoords."
        isa(neighbors,Int) && (neighbors=Dict(i=>minlen for (i,minlen) in enumerate(minimumlengths(rcoords,vectors,neighbors,coordination=coordination))))
        vectors=convert(Vector{SVector{N,Float}},vectors)
        recipls=convert(Vector{SVector{N,Float}},reciprocals(vectors))
        new{N,pids|>eltype}(name,pids,rcoords,icoords,vectors,recipls,neighbors)
    end
end
function Lattice(   name::String,
                    points::AbstractVector{<:Point};
                    vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{points|>eltype|>dimension,Float}}(),
                    neighbors::Union{Dict{Int,<:Real},Int}=1,
                    coordination::Int=8
                    )
    pids=Vector{points|>eltype|>pidtype}(undef,points|>length)
    rcoords=zeros(Float,points|>eltype|>dimension,points|>length)
    icoords=zeros(Float,points|>eltype|>dimension,points|>length)
    for i=1:length(points)
        pids[i]=points[i].pid
        rcoords[:,i]=points[i].rcoord
        icoords[:,i]=points[i].icoord
    end
    Lattice{points|>eltype|>dimension}(name,pids,rcoords,icoords,vectors,neighbors,coordination=coordination)
end
function Lattice(   name::String,
                    sublattices::AbstractVector{<:AbstractLattice};
                    vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                    neighbors::Union{Dict{Int,<:Real},Int}=1,
                    coordination::Int=8
                    )
    if length(sublattices)>1
        @assert mapreduce(sl->dimension(sl)==dimension(sublattices[1]),&,sublattices) "Lattice error: all sublattices should have the same dimension."
        @assert mapreduce(sl->keytype(sl)===keytype(sublattices[1]),&,sublattices) "Lattice error: all sublattices should have the same keytype."
    end
    len=mapreduce(length,+,sublattices)
    pids=Vector{sublattices|>eltype|>keytype}(undef,len)
    rcoords=zeros(Float,sublattices|>eltype|>dimension,len)
    icoords=zeros(Float,sublattices|>eltype|>dimension,len)
    pos=1
    for sublattice in sublattices
        pids[pos:pos+length(sublattice)-1]=sublattice.pids
        rcoords[:,pos:pos+length(sublattice)-1]=sublattice.rcoords
        icoords[:,pos:pos+length(sublattice)-1]=sublattice.icoords
        pos+=length(sublattice)
    end
    Lattice{sublattices|>eltype|>dimension}(name,pids,rcoords,icoords,vectors,neighbors,coordination=coordination)
end

"""
    SuperLattice(   name::String,
                    sublattices::AbstractVector{<:AbstractLattice};
                    vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                    neighbors::Dict{Int,<:Real}=Dict{Int,Float}()
                    )

SuperLattice that is composed of serveral sublattices.
"""
struct SuperLattice{L<:AbstractLattice,N,P<:PID} <: AbstractLattice{N,P}
    sublattices::Vector{L}
    name::String
    pids::Vector{P}
    rcoords::Matrix{Float}
    icoords::Matrix{Float}
    vectors::Vector{SVector{N,Float}}
    reciprocals::Vector{SVector{N,Float}}
    neighbors::Dict{Int,Float}
    function SuperLattice(  name::String,
                            sublattices::AbstractVector{<:AbstractLattice};
                            vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                            neighbors::Dict{Int,<:Real}=Dict{Int,Float}()
                            )
        if length(sublattices)>1
            @assert mapreduce(sl->dimension(sl)==dimension(sublattices[1]),&,sublattices) "SuperLattice error: all sublattices should have the same dimension."
            @assert mapreduce(sl->keytype(sl)===keytype(sublattices[1]),&,sublattices) "SuperLattice error: all sublattices should have the same keytype."
        end
        @assert all(length(sublattice.vectors)==0 for sublattice in sublattices) "SuperLattice error: all sublattices should assume open boundaries."
        len=mapreduce(length,+,sublattices)
        pids=Vector{sublattices|>eltype|>keytype}(undef,len)
        rcoords=zeros(Float,sublattices|>eltype|>dimension,len)
        icoords=zeros(Float,sublattices|>eltype|>dimension,len)
        pos=1
        for sublattice in sublattices
            pids[pos:pos+length(sublattice)-1]=sublattice.pids
            rcoords[:,pos:pos+length(sublattice)-1]=sublattice.rcoords
            icoords[:,pos:pos+length(sublattice)-1]=sublattice.icoords
            pos+=length(sublattice)
        end
        vectors=convert(Vector{SVector{sublattices|>eltype|>dimension,Float}},vectors)
        recipls=convert(Vector{SVector{sublattices|>eltype|>dimension,Float}},reciprocals(vectors))
        new{sublattices|>eltype,sublattices|>eltype|>dimension,sublattices|>eltype|>keytype}(sublattices,name,pids,rcoords,icoords,vectors,recipls,neighbors)
    end
end

"""
    latticetype(sl::SuperLattice)
    latticetype(::Type{<:SuperLattice{L}}) where L<:AbstractLattice

Get the sublattice type of a superlattice.
"""
latticetype(sl::SuperLattice)=sl|>typeof|>latticetype
latticetype(::Type{<:SuperLattice{L}}) where L<:AbstractLattice=L

struct IntraBonds <: LatticeBonds{2} end
struct InterBonds <: LatticeBonds{2} end
"""
    intrabonds

Indicate that bonds intra the sublattices are inquired.
!!! note
    These bonds do not contain those accorss the periodic boundaries.
"""
const intrabonds=IntraBonds()
"""
    interbonds

Indicate that bonds inter the sublattices are inquired.
!!! note
    These bonds do not contain those accorss the periodic boundaries.
"""
const interbonds=InterBonds()

"""
    latticebondsstructure(::Type{<:SuperLattice}) -> SimpleTree{LatticeBonds,Nothing}

The tree structure of the lattice bonds.
"""
@generated function latticebondsstructure(::Type{<:SuperLattice})
    structure=SimpleTree{LatticeBonds,Nothing}()
    push!(structure,allbonds,nothing)
    push!(structure,allbonds,zerothbonds,nothing)
    push!(structure,allbonds,insidebonds,nothing)
    push!(structure,allbonds,acrossbonds,nothing)
    push!(structure,insidebonds,intrabonds,nothing)
    push!(structure,insidebonds,interbonds,nothing)
    return structure
end

"""
    bonds!(bonds::Vector,lattice::SuperLattice,::Val{intrabonds}) -> Vector
    bonds!(bonds::Vector,lattice::SuperLattice,::Val{interbonds}) -> Vector

Generate the required bonds of a superlattice and append them to the input bonds.
"""
function bonds!(bonds::Vector,lattice::SuperLattice,::Val{intrabonds})
    for sublattice in lattice.sublattices
        bonds!(bonds,sublattice,insidebonds)
    end
    return bonds
end
function bonds!(bonds::Vector,lattice::SuperLattice,::Val{interbonds})
    dim=lattice|>dimension|>Val
    for (i,j) in Combinations{2}(1:length(lattice.sublattices))
        sub1,sub2=lattice.sublattices[i],lattice.sublattices[j]
        for (neighbor,sindex,eindex) in interlinks(sub1.rcoords,sub2.rcoords,lattice.neighbors)
            spoint=Point(sub1.pids[sindex],latticestaticcoords(sub1.rcoords,sindex,dim),latticestaticcoords(sub1.icoords,sindex,dim))
            epoint=Point(sub2.pids[eindex],latticestaticcoords(sub2.rcoords,eindex,dim),latticestaticcoords(sub2.icoords,eindex,dim))
            push!(bonds,Bond(neighbor,spoint,epoint))
        end
    end
    return bonds
end

"""
    Cylinder{P}(name::String,
                block::AbstractMatrix{<:Real},
                translation::SVector{N,<:Real};
                vector::Union{AbstractVector{<:Real},Nothing}=nothing,
                neighbors::Union{Dict{Int,<:Real},Int}=1
                ) where {P<:PID,N}

Cylinder of 1d and quasi 2d lattices.
"""
mutable struct Cylinder{P<:PID,N} <: AbstractLattice{N,P}
    block::Matrix{Float}
    translation::SVector{N,Float}
    name::String
    pids::Vector{P}
    rcoords::Matrix{Float}
    icoords::Matrix{Float}
    vectors::Vector{SVector{N,Float}}
    reciprocals::Vector{SVector{N,Float}}
    neighbors::Dict{Int,Float}
    function Cylinder{P}(   name::String,
                            block::AbstractMatrix{<:Real},
                            translation::SVector{N,<:Real};
                            vector::Union{AbstractVector{<:Real},Nothing}=nothing,
                            neighbors::Union{Dict{Int,<:Real},Int}=1
                            ) where {P<:PID,N}
        pids=Vector{P}[]
        rcoords=zeros(Float,N,0)
        icoords=zeros(Float,N,0)
        vectors=vector===nothing ? SVector{N,Float}[] : [convert(SVector{N,Float},vector)]
        recipls=convert(Vector{SVector{N,Float}},reciprocals(vectors))
        isa(neighbors,Int) && (neighbors=Dict{Int,Float}(i=>Inf for i=1:neighbors))
        new{P,N}(block,translation,name,pids,rcoords,icoords,vectors,recipls,neighbors)
    end
end

"""
    insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=9) where S -> Cylinder

Insert a couple of blocks into a cylinder.

The position of the cut of the cylinder is specified by the keyword argument `cut`, which is the center of the cylinder by default. All pids corresponding to a same newly inserted block share the same scope, which is specified by the parameter `ps`. Optionally, the scopes of the old pids in the cylinder can be replaced if the parameter `scopes` is assigned other than `nothing`. Note the length of `ps` is equal to the number of newly inserted blocks, while that of `scopes` should be equal to the old length of the cylinder.
"""
function Base.insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=8) where S
    @assert S<:fieldtype(cylinder|>keytype,:scope) "insert! error: dismatched type of input scopes and old scopes."
    @assert 1<=cut<=length(cylinder)+1 "insert! error: wrong cut($cut), which should be in [1,$(length(cylinder)+1)]."
    @assert length(cylinder)%size(cylinder.block,2)==0 "insert! error: wrong cylinder length."
    scopes===nothing || @assert length(scopes)==length(cylinder) "insert! error: dismatched length of input scopes and cylinder."
    blcklen=size(cylinder.block,2)
    blcknum=length(cylinder)÷blcklen+length(ps)
    rsltlen=blcklen*blcknum
    pids=keytype(cylinder)[]
    rcoords=zeros(Float,cylinder|>dimension,rsltlen)
    icoords=zeros(Float,cylinder|>dimension,rsltlen)
    for i=1:cut-1
        push!(pids,scopes===nothing ? cylinder.pids[i] : replace(cylinder.pids[i],scope=scopes[i]))
    end
    for (i,p) in enumerate(ps)
        for j=1:blcklen
            push!(pids,PID(p,cut==length(cylinder)+1 ? j : (cylinder.pids[cut].site+j-2)%blcklen+1))
        end
    end
    for i=cut:length(cylinder)
        push!(pids,scopes===nothing ? cylinder.pids[i] : replace(cylinder.pids[i],scope=scopes[i]))
    end
    for i=1:blcknum
        trannum=i-(blcknum+1)/2
        for j=1:blcklen
            for k=1:dimension(cylinder)
                rcoords[k,(i-1)*blcklen+j]=cylinder.block[k,j]+trannum*cylinder.translation[k]
            end
        end
    end
    cylinder.pids=pids
    cylinder.rcoords=rcoords
    cylinder.icoords=icoords
    if length(cylinder.neighbors)>0 && mapreduce(isequal(Inf),|,values(cylinder.neighbors))
        minlens=minimumlengths(rcoords,cylinder.vectors,cylinder|>nneighbor,coordination=coordination)
        cylinder.neighbors=Dict{Int,Float}(i=>len for (i,len) in enumerate(minlens))
    end
    cylinder
end

"""
    (cylinder::Cylinder)(scopes::Any...;coordination::Int=8) -> Lattice

Construct a lattice from a cylinder with the assigned scopes.
"""
function (cylinder::Cylinder)(scopes::Any...;coordination::Int=8)
    @assert fieldtype(cylinder|>keytype,:scope)===scopes|>eltype "cylinder call error: wrong scope type."
    name=cylinder.name*string(length(scopes))
    pids=keytype(cylinder)[PID(scope,i) for scope in scopes for i=1:size(cylinder.block,2)]
    rcoords=tile(cylinder.block,[cylinder.translation],Tuple((i,) for i=-(length(scopes)-1)/2:(length(scopes)-1)/2))
    icoords=zero(rcoords)
    vectors=cylinder.vectors
    neighbors=cylinder.neighbors
    if length(neighbors)>0 && mapreduce(isequal(Inf),|,values(neighbors))
        minlens=minimumlengths(rcoords,vectors,cylinder|>nneighbor,coordination=coordination)
        neighbors=Dict{Int,Float}(i=>len for (i,len) in enumerate(minlens))
    end
    Lattice{cylinder|>dimension}(name,pids,rcoords,icoords,vectors,neighbors)
end

"""
    Bonds{T,L}(bonds::Tuple{Vararg{Vector{<:AbstractBond}}}) where {T,L<:AbstractLattice}
    Bonds(lattice::AbstractLattice,types::LatticeBonds...)

A set of lattice bonds.

`Bonds` itself is an `AbstractVector` of `AbstractBond`. The need for such a struct is to ensure the type stability during the iteration over a set of different concrete bonds. Although the default `iterate` function does not achieve this goal, users can get it with the generated function trick. Besides, it provides a high level of management of different categories of bonds based on the `LatticeBonds` system.
"""
struct Bonds{T,L<:AbstractLattice,BS<:Tuple{Vararg{Vector{<:AbstractBond}}},B<:AbstractBond} <: AbstractVector{B}
    bonds::BS
    function Bonds{T,L}(bonds::Tuple{Vararg{Vector{<:AbstractBond}}}) where {T,L<:AbstractLattice}
        @assert isa(T,Tuple{Vararg{LatticeBonds}}) && length(T)==length(bonds) "Bonds error: dismatched input types and bonds."
        new{T,L,typeof(bonds),mapreduce(eltype,typejoin,bonds)}(bonds)
    end
end
Bonds(lattice::AbstractLattice)=Bonds(lattice,allbonds)
Bonds(lattice::AbstractLattice,types::LatticeBonds...)=Bonds(lattice,Val(types))
@generated function Bonds(lattice::AbstractLattice,::Val{VS}) where VS
    types,bonds=[],[]
    for V in expand(lattice,Val(VS))
        push!(types,V)
        push!(bonds,:(bonds(lattice,Val($V))))
    end
    types=Tuple(types)
    bonds=Expr(:tuple,bonds...)
    return :(Bonds{$types,typeof(lattice)}($bonds))
end

"""
    ==(bonds1::Bonds,bonds2::Bonds) -> Bool
    isequal(bonds1::Bonds,bonds2::Bonds) -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(bonds1::Bonds,bonds2::Bonds) = ==(efficientoperations,bonds1,bonds2)
Base.isequal(bonds1::Bonds,bonds2::Bonds)=isequal(efficientoperations,bonds1,bonds2)

"""
    size(bs::Bonds) -> Tuple{Int}

Get the size of the set of lattice bonds.
"""
Base.size(bs::Bonds)=(mapreduce(length,+,bs.bonds),)

"""
    summary(io::IO,bs::Bonds)

Print the brief description of a set of lattice bonds to an io.
"""
Base.summary(io::IO,bs::Bonds)=@printf io "%s-element Bonds" length(bs)

"""
    bondtypes(bs::Bonds) -> Tuple{Vararg{LatticeBonds}}
    bondtypes(::Type{<:Bonds{T}}) where T -> Tuple{Vararg{LatticeBonds}}

Get the bondtypes of a set of lattice bonds.
"""
bondtypes(bs::Bonds)=bs|>typeof|>bondtypes
bondtypes(::Type{<:Bonds{T}}) where T=T

"""
    latticetype(bs::Bonds)
    latticetype(::Type{<:Bonds{T,L} where T}) where L<:AbstractLattice
"""
latticetype(bs::Bonds)=bs|>typeof|>latticetype
latticetype(::Type{<:Bonds{T,L} where T}) where L<:AbstractLattice=L

"""
    rank(bs::Bonds) -> Int
    rank(::Type{<:Bonds{T}}) where T -> Int

Get the rank of a set of lattice bonds.
"""
rank(bs::Bonds)=bs|>typeof|>rank
rank(::Type{<:Bonds{T}}) where T=length(T)

"""
    iterate(bs::Bonds,state=(1,0))

Iterate over the lattice bonds in the set.
"""
function Base.iterate(bs::Bonds,state=(1,0))
    s1,s2=state[1],state[2]
    while s1<=length(bs.bonds)
        s1,s2=s2+1>length(bs.bonds[s1]) ? (s1+1,1) : (s1,s2+1)
        (s1>length(bs.bonds) || s2<=length(bs.bonds[s1])) && break
    end
    s1>length(bs.bonds) && return nothing
    return bs.bonds[s1][s2],(s1,s2)
end

"""
    getindex(bs::Bonds,i::Int) -> eltype(bs)

Get the ith bond in the set.
"""
function Base.getindex(bs::Bonds,i::Int)
    r,k=1,i
    while r<=rank(bs) && k>length(bs.bonds[r])
        k=k-length(bs.bonds[r])
        r=r+1
    end
    r>rank(bs) && error("getindex error: attempt to access $(length(bs))-element Bonds at index [$k].")
    return bs.bonds[r][k]
end

"""
    filter(lbs::LatticeBonds,bs::Bonds,choice::Union{Val{:include},Val{:exclude}}=Val(:include)) -> Bonds
    filter(lbs::Tuple{Vararg{LatticeBonds}},bs::Bonds,choice::Union{Val{:include},Val{:exclude}}=Val(:include)) -> Bonds

Get a subset of a set of lattice bonds.

When `choice=Val(:include)`, the lattice bonds indicated by `lbs` will be selected;
When `choice=Val(:exclude)`, the lattice bonds not indicated by `lbs` will be selected.
"""
Base.filter(lbs::LatticeBonds,bs::Bonds,choice::Union{Val{:include},Val{:exclude}}=Val(:include))=bondfilter(Val((lbs,)),bs,choice)
Base.filter(lbs::Tuple{Vararg{LatticeBonds}},bs::Bonds,choice::Union{Val{:include},Val{:exclude}}=Val(:include))=bondfilter(Val(lbs),bs,choice)
@generated function bondfilter(::Val{VS},bs::Bonds,::Val{:include}) where VS
    types,bonds=[],[]
    for V in expand(latticetype(bs),Val(VS))
        index=findfirst(isequal(V),bondtypes(bs))
        if isa(index,Int)
            push!(types,V)
            push!(bonds,:(bs.bonds[$index]))
        end
    end
    types=Tuple(types)
    bonds=Expr(:tuple,bonds...)
    return :(Bonds{$types,latticetype(bs)}($bonds))
end
@generated function bondfilter(::Val{VS},bs::Bonds,::Val{:exclude}) where VS
    types,bonds=[],[]
    excludes=expand(latticetype(bs),Val(VS))
    for (i,V) in enumerate(bondtypes(bs))
        if V ∉ excludes
            push!(types,V)
            push!(bonds,:(bs.bonds[$i]))
        end
    end
    types=Tuple(types)
    bonds=Expr(:tuple,bonds...)
    return :(Bonds{$types,latticetype(bs)}($bonds))
end

"""
    filter(select::Function,bs::Bonds) -> Bonds

Get a filtered set of bonds by a select function.
"""
@generated function Base.filter(select::Function,bs::Bonds)
    bonds=Expr(:tuple,[:(filter(select,bs.bonds[$i])) for i=1:rank(bs)]...)
    return :(Bonds{bondtypes(bs),latticetype(bs)}($bonds))
end

"""
    empty!(bs::Bonds) -> Bonds

Empty a set of lattice bonds.
"""
@generated function Base.empty!(bs::Bonds)
    exprs=[:(empty!(bs.bonds[$i])) for i=1:rank(bs)]
    push!(exprs,:(return bs))
    return Expr(:block,exprs...)
end

"""
    empty(bs::Bonds) -> Bonds

Get an empty copy of a set of lattice bonds.
"""
@generated function Base.empty(bs::Bonds)
    bonds=Expr(:tuple,[:(empty(bs.bonds[$i])) for i=1:rank(bs)]...)
    return :(Bonds{bondtypes(bs),latticetype(bs)}($bonds))
end

"""
    reset!(bs::Bonds,lattice::AbstractLattice) -> Bonds

Reset a set of lattice bonds by a new lattice.
"""
@generated function reset!(bs::Bonds,lattice::AbstractLattice)
    exprs=[]
    push!(exprs,quote
        @assert latticetype(bs)>:typeof(lattice) "reset! error: dismatched bonds and lattice."
        empty!(bs)
    end)
    for i=1:rank(bs)
        push!(exprs,:(bonds!(bs.bonds[$i],lattice,bondtypes(bs)[$i])))
    end
    push!(exprs,:(return bs))
    return Expr(:block,exprs...)
end

end #module
