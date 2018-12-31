module Spatial

using LinearAlgebra: norm,dot,cross,inv
using StaticArrays: SVector,SMatrix
using Printf: @printf,@sprintf
using NearestNeighbors: KDTree,knn,inrange
using Base.Iterators: flatten,product
using ...Utilities: atol,rtol,Float
using ...Utilities.TypeTrait: efficientoperations
using ...Utilities.Combinatorics: Combinations
using ...Utilities.NamedVector: AbstractNamedVector
using ...Utilities.Factory: Inference,TypeFactory,FunctionFactory,Argument,MixEscaped,Escaped,UnEscaped
using ...Utilities.Factory: addparams!,addfields!,addwhereparams!,addargs!,addkwargs!,extendbody!,addconstructors!

import ...Utilities.Interface: rank,dimension

export rank,dimension
export distance,azimuthd,azimuth,polard,polar,volume
export isparallel,isonline,isintratriangle,issubordinate
export reciprocals,translate,rotate,tile,minimumlengths
export Link,intralinks,interlinks
export PID,AbstractBond,Point,Bond,pidtype,rcoord,icoord,isintracell
export AbstractLattice,nneighbor,bonds
export AbstractLatticeIndex,RCoordIndex,ICoordIndex,PointIndex
export zerothbonds,insidebonds,acrossbonds,intrabonds,interbonds
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
    result=acosd(v[1]/(norm(length(v)==3 ? @views(v[1:end-1]) : v)))
    length(v)>1 && v[2]<0 && (result=360-result)
    return result
end

"""
    azimuth(v::AbstractVector{<:Real}) -> Float

Get the azimuth angle in radians of a vector.
"""
function azimuth(v::AbstractVector{<:Real})
    @assert length(v) ∈ (1,2,3) "azimuth error: wrong dimensioned input vector."
    result=acos(v[1]/(norm(length(v)==3 ? @views(v[1:end-1]) : v)))
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
        result=dot(v1,cross(v2,v3))
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
!!! notes
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
    ndim=length(p)
    a,b=zeros(Float,3,3),zeros(Float,3)
    a[1:ndim,1]=p2-p1
    a[1:ndim,2]=p3-p1
    a[:,3]=cross(a[:,1],a[:,2])
    b[1:ndim]=p-p1
    x=inv(a)*b
    @assert isapprox(x[3],0.0,atol=atol,rtol=rtol) "isintratriangle error: internal error..."
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
    @assert all(length(vector)==length(rcoord) for vector in vectors) "issubordinate error: shape dismatch of input rcoord and vectors."
    ndim,nv=length(rcoord),length(vectors)
    a,b=zeros(Float,3,3),zeros(Float,3)
    for i=1:nv
        a[1:ndim,i]=vectors[i]
    end
    if nv==1 || nv==2
        if nv==1
            for i=1:3
                a[:,2].=0.0
                a[i,2]=pi
                isparallel(@views(a[:,1]),@views(a[:,2]),atol=atol,rtol=rtol)==0 && break
            end
        end
        a[:,3]=cross(@views(a[:,1]),@views(a[:,2]))
    end
    b[1:ndim]=rcoord
    x=inv(a)*b
    return isapprox(norm(round.(x)-x),0.0,atol=atol,rtol=rtol)
end

"""
    reciprocals(vectors::AbstractVector{AbstractVector{<:Real}}) -> Vector{Vector{Float}}

Get the reciprocals dual to the input vectors.
"""
function reciprocals(vectors::AbstractVector{<:AbstractVector{<:Real}})
    @assert length(vectors)<4 "reciprocals error: the number of input vectors should not be greater than 3."
    @assert all(length(vector) in (1,2,3) for vector in vectors) "reciprocals error: all input vectors must be 1, 2 or 3 dimensional."
    result=Vector{Float}[]
    if length(vectors)==1
        push!(result,2pi*vectors[1]/sum(vectors[1].^2))
    elseif length(vectors) in (2,3)
        ndim=length(vectors[1])
        m=zeros(Float,3,3)
        m[1:ndim,1]=vectors[1]
        m[1:ndim,2]=vectors[2]
        m[:,3]=length(vectors)==2 ? cross(@views(m[:,1]),@views(m[:,2])) : vectors[3]
        m=inv(m)
        push!(result,2pi*m[1,1:ndim])
        push!(result,2pi*m[2,1:ndim])
        length(vectors)==3 && push!(result,2*pi*m[3,1:ndim])
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
!!! notes
    1. The result is given by the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula).
    2. Only 2 and 3 dimensional vectors can be rotated.
    3. When the input vectors are 2 dimensional, both the polar and azimuth of the axis must be 0.
"""
function rotate(cluster::AbstractMatrix{<:Real},angle::Real;axis::Tuple{Union{AbstractVector{<:Real},Nothing},Tuple{<:Real,<:Real}}=(nothing,(0,0)))
    @assert size(cluster,1) in (2,3) "rotate error: only 2 and 3 dimensional vectors can be rotated."
    center,theta,phi=(axis[1]===nothing ? zeros(size(cluster,1)) : axis[1]),axis[2][1],axis[2][2]
    @assert length(center)==size(cluster,1) "rotate error: dismatched shape of the input cluster and the point on axis."
    length(center)==2 && @assert isapprox(theta,0,atol=atol) && isapprox(phi,0,atol=atol) "rotate error: both the polar and azimuth of the axis for 2 dimensional vectors must be 0."
    cache,cosθ,sinθ=zeros(3),cos(angle),sin(angle)
    k=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
    result=zeros(size(cluster))
    for i=1:size(cluster,2)
        for j=1:size(cluster,1)
            cache[j]=cluster[j,i]-center[j]
        end
        inner,outer=dot(k,cache),cross(k,cache)
        for j=1:size(cluster,1)
            result[j,i]=cache[j]*cosθ+outer[j]*sinθ+k[j]*inner*(1-cosθ)+center[j]
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
    N==0 || @assert length(vectors)==N "tile error: dismatched shape of input vectors and translations."
    supercluster=zeros(size(cluster,1),size(cluster,2)*length(translations))
    for (i,translation) in enumerate(N==0 ? (0,) : translations)
        disp=N==0 ? zeros(size(cluster,1)) : sum(vectors.*translation)
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
            if any(translation.!=0)
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
        any(result.==Inf) && @warn "minimumlengths warning: Inf remained in the result and larger(>$coordination) coordination or smaller(<$nneighbor) nneighbor may be needed."
    end
    return result
end

"""
    Link(neighbor::Int,sindex::Int,eindex::Int,disp::AbstractVector{<:Real})

A link in a lattice.
"""
struct Link{N}
    neighbor::Int
    sindex::Int
    eindex::Int
    disp::SVector{N,Float}
end
Link(neighbor::Int,sindex::Int,eindex::Int,disp::AbstractVector{<:Real})=Link(neighbor,sindex,eindex,convert(SVector{disp|>length,Float},disp))

"""
    show(io::IO,link::Link)

Show a link.
"""
Base.show(io::IO,link::Link)=@printf io "Link(%s,%s,%s,[%s])" link.neighbor link.sindex link.eindex join(string.(link.disp),",")

"""
    ==(l1::Link,l2::Link) -> Bool
    isequal(l1::Link,l2::Link) -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(l1::Link,l2::Link) = ==(efficientoperations,l1,l2)
Base.isequal(l1::Link,l2::Link)=isequal(efficientoperations,l1,l2)

"""
    dimension(link::Link) -> Int
    dimension(::Type{<:Link{N}}) where N -> Int

Get the space dimension of a link.
"""
dimension(link::Link)=link|>typeof|>dimension
dimension(::Type{<:Link{N}}) where N=N

"""
    intralinks( cluster::AbstractMatrix{<:Real},
                vectors::AbstractVector{<:AbstractVector{<:Real}},
                neighbors::Dict{Int,Float},
                maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))
                ) where N -> Vector{Link}

Use kdtree to get the intracluster nearest neighbors.

As is similar to [`minimumlengths`](@ref), when `vectors` is nonempty, the cluster assumes periodic boundaries. `neighbors` provides the map between the bond length and the order of nearest neighbors. Note only those with the lengths present in `neighbors` will be included in the result. `maxtranslations` determines the maximum number of translations along those directions specified by `vectors` when the tiled supercluster is construted (See [`minimumlengths`](@ref) for the explanation of the method for periodic lattices).
"""
function intralinks(    cluster::AbstractMatrix{<:Real},
                        vectors::AbstractVector{<:AbstractVector{<:Real}},
                        neighbors::Dict{Int,Float},
                        maxtranslations::NTuple{N,Int}=ntuple(i->length(neighbors),length(vectors))
                        ) where N
    @assert length(vectors)==N "intralinks error: dismatched shape of input vectors and maxtranslations."
    result=Link{size(cluster,1)}[]
    translations=reshape(product((-nnb:nnb for nnb in maxtranslations)...)|>collect,:)
    for translation in translations
        if any(translation.!=0)
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
                dist=norm(@views(supercluster[:,j]-cluster[:,i]))
                for (nb,len) in neighbors
                    isapprox(len,dist,atol=atol) && (push!(result,Link(nb,i,(j-1)%size(cluster,2)+1,disps[:,j]));break)
                end
            end
        end
    end
    return result
end

"""
    interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float}) -> Vector{Link}

Use kdtree to get the intercluster nearest neighbors.
"""
function interlinks(cluster1::AbstractMatrix{<:Real},cluster2::AbstractMatrix{<:Real},neighbors::Dict{Int,Float})
    @assert size(cluster1,1)==size(cluster1,1) "interlinks error: dismatched space dimension of input clusters."
    result=Link{size(cluster1,1)}[]
    for (i,indices) in enumerate(inrange(KDTree(cluster2),cluster1,max(values(neighbors)...)+atol,true))
        for j in indices
            dist=norm(@views(cluster2[:,j]-cluster1[:,i]))
            for (nb,len) in neighbors
                isapprox(len,dist,atol=atol) && (push!(result,Link(nb,i,j,zero(SVector{size(cluster1,1),Float})));break)
            end
        end
    end
    return result
end

"""
    PID(scope,site::Int)
    PID(;scope="tz",site::Int=1)

The id of a point.
"""
struct PID{S} <: AbstractNamedVector
    scope::S
    site::Int
end
PID(;scope="tz",site::Int=1)=PID(scope,site)
Base.fieldnames(pid::PID)=(:scope,:site)

"""
    AbstractBond{R,P<:PID,N}

Abstract bond.
"""
abstract type AbstractBond{R,P<:PID,N} end

"""
    ==(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool
    isequal(b1::AbstractBond{R},b1::AbstractBond{R}) where R -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(b1::AbstractBond{R},b2::AbstractBond{R}) where R = ==(efficientoperations,b1,b2)
Base.isequal(b1::AbstractBond{R},b2::AbstractBond{R}) where R=isequal(efficientoperations,b1,b2)

"""
    rank(b::AbstractBond) -> Int
    rank(::Type{<:AbstractBond{R,P,N}}) where {R,P,N} -> Int

Get the rank of a concrete bond.
"""
rank(b::AbstractBond)=b|>typeof|>rank
rank(::Type{<:AbstractBond{R,P,N}}) where {R,P,N}=R

"""
    pidtype(b::AbstractBond)
    pidtype(::Type{<:AbstractBond{R,P,N}}) where {R,P,N}

Get the pid type of a concrete bond.
"""
pidtype(b::AbstractBond)=b|>typeof|>pidtype
pidtype(::Type{<:AbstractBond{R,P,N}}) where {R,P,N}=P

"""
    dimension(b::AbstractBond) -> Int
    dimension(::Type{<:AbstractBond{R,P,N}}) where {R,P,N} -> Int

Get the space dimension of a concrete bond.
"""
dimension(b::AbstractBond)=b|>typeof|>dimension
dimension(::Type{<:AbstractBond{R,P,N}}) where {R,P,N}=N

"""
    Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N
    Point(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N
    Point(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))

Labeled point.
"""
struct Point{P<:PID,N} <: AbstractBond{1,P,N}
    pid::P
    rcoord::SVector{N,Float}
    icoord::SVector{N,Float}
    Point(pid::PID,rcoord::SVector{N,<:Real},icoord::SVector{N,<:Real}) where N=new{pid|>typeof,N}(pid,rcoord,icoord)
end
Point(pid::PID,rcoord::NTuple{N,<:Real},icoord::NTuple{N,<:Real}=ntuple(i->0.0,N)) where N=Point(pid,convert(SVector{N,Float},rcoord),convert(SVector{N,Float},icoord))
function Point(pid::PID,rcoord::AbstractVector{<:Real},icoord::AbstractVector{<:Real}=zero(SVector{length(rcoord),Float}))
    @assert length(rcoord)==length(icoord) "Point error: dismatched length of input rcoord and icoord."
    Point(pid,convert(SVector{length(rcoord),Float},rcoord),convert(SVector{length(icoord),Float},icoord))
end

"""
    show(io::IO,p::Point)

Show a labeled point.
"""
Base.show(io::IO,p::Point)=@printf io "Point(%s,[%s],[%s])" p.pid join(string.(p.rcoord),",") join(string.(p.icoord),",")

"""
    Bond(neighbor::Int,spoint::Point,epoint::Point)

A bond in a lattice.
"""
struct Bond{P<:PID,N} <: AbstractBond{2,P,N}
    neighbor::Int
    spoint::Point{P,N}
    epoint::Point{P,N}
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
abstract type AbstractLattice{P<:PID,N} end

"""
    show(io::IO,lattice::AbstractLattice)

Show a lattice.
"""
Base.show(io::IO,lattice::AbstractLattice)=@printf io "%s(%s)" lattice|>typeof|>nameof lattice.name

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
    dimension(lattice::AbstractLattice) -> Int
    dimension(::Type{<:AbstractLattice{P,N}}) where {P,N}-> Int

Get the space dimension of the lattice.
"""
dimension(lattice::AbstractLattice)=lattice|>typeof|>dimension
dimension(::Type{<:AbstractLattice{P,N}}) where {P,N}=N

"""
    keytype(lattice::AbstractLattice)
    keytype(::Type{<:AbstractLattice{P,N}}) where {P,N}

Get the pid type of the lattice.
"""
Base.keytype(lattice::AbstractLattice)=lattice|>typeof|>keytype
Base.keytype(::Type{<:AbstractLattice{P,N}}) where {P,N}=P

"""
    valtype(lattice::AbstractLattice)
    valtype(::Type{<:AbstractLattice{P,N}}) where {P,N}

Get the point type of the lattice.
"""
Base.valtype(lattice::AbstractLattice)=lattice|>typeof|>valtype
Base.valtype(::Type{<:AbstractLattice{P,N}}) where {P,N}=Point{P,N}

"""
    AbstractLatticeIndex{I<:Union{<:PID,Int}}

Abstract index type for a lattice.
"""
abstract type AbstractLatticeIndex{I<:Union{<:PID,Int}} end

"""
    RCoordIndex(index::Union{<:PID,Int})

Index for getting a rcoord of a lattice.
"""
struct RCoordIndex{I} <: AbstractLatticeIndex{I} index::I end

"""
    ICoordIndex(index::Union{<:PID,Int})

Index for getting an icoord of a lattice.
"""
struct ICoordIndex{I} <: AbstractLatticeIndex{I} index::I end

"""
    PointIndex(index::Union{<:PID,Int})

Index for getting a point of a lattice.
"""
struct PointIndex{I} <: AbstractLatticeIndex{I} index::I end

"""
    getindex(lattice::AbstractLattice,i::RCoordIndex) -> SVector
    getindex(lattice::AbstractLattice,i::ICoordIndex) -> SVector
    getindex(lattice::AbstractLattice,i::PointIndex) -> Point

Get a rcoord, an icoord or a point of a lattice according to the type of the input index.
"""
Base.getindex(lattice::AbstractLattice,i::RCoordIndex{Int})=SVector{lattice|>dimension}(@views(lattice.rcoords[:,i.index]))
Base.getindex(lattice::AbstractLattice,i::ICoordIndex{Int})=SVector{lattice|>dimension}(@views(lattice.icoords[:,i.index]))
Base.getindex(lattice::AbstractLattice,i::PointIndex{Int})=@views Point(lattice.pids[i.index],lattice.rcoords[:,i.index],lattice.icoords[:,i.index])
Base.getindex(lattice::AbstractLattice,i::RCoordIndex{<:PID})=lattice[RCoordIndex(findfirst(isequal(i.index),lattice.pids))]
Base.getindex(lattice::AbstractLattice,i::ICoordIndex{<:PID})=lattice[ICoordIndex(findfirst(isequal(i.index),lattice.pids))]
Base.getindex(lattice::AbstractLattice,i::PointIndex{<:PID})=lattice[PointIndex(findfirst(isequal(i.index),lattice.pids))]

"""
    nneighbor(lattice::AbstractLattice) -> Int

Get the highest order of nearest neighbors.
"""
nneighbor(lattice::AbstractLattice)=max(keys(lattice.neighbors)...)

abstract type LatticeBonds end
struct ZerothBonds <: LatticeBonds end
struct InsideBonds <: LatticeBonds end
struct AcrossBonds <: LatticeBonds end
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
    bonds(lattice::AbstractLattice) -> Vector{AbstractBond}
    bonds(lattice::AbstractLattice,::ZerothBonds) -> Vector{Point}
    bonds(lattice::AbstractLattice,::InsideBonds) -> Vector{Bond}
    bonds(lattice::AbstractLattice,::AcrossBonds) -> Vector{Bond}

Get the bonds of a lattice.
"""
function bonds(lattice::AbstractLattice)
    result=convert(Vector{AbstractBond},bonds(lattice,zerothbonds))
    for link in intralinks(lattice.rcoords,lattice.vectors,lattice.neighbors)
        @views spoint=Point(lattice.pids[link.sindex],lattice.rcoords[:,link.sindex],lattice.icoords[:,link.sindex])
        @views epoint=Point(lattice.pids[link.eindex],lattice.rcoords[:,link.eindex]+link.disp,lattice.icoords[:,link.eindex]+link.disp)
        push!(result,Bond(link.neighbor,spoint,epoint))
    end
    result
end
function bonds(lattice::AbstractLattice,::ZerothBonds)
    return [Point(lattice.pids[i],@views(lattice.rcoords[:,i]),@views(lattice.icoords[:,i])) for i=1:length(lattice)]
end
function bonds(lattice::AbstractLattice,::InsideBonds)
    result=Bond{lattice|>keytype,lattice|>dimension}[]
    for link in interlinks(lattice.rcoords,lattice.rcoords,lattice.neighbors)
        if link.sindex<=link.eindex
            @views spoint=Point(lattice.pids[link.sindex],lattice.rcoords[:,link.sindex],lattice.icoords[:,link.sindex])
            @views epoint=Point(lattice.pids[link.eindex],lattice.rcoords[:,link.eindex]+link.disp,lattice.icoords[:,link.eindex]+link.disp)
            push!(result,Bond(link.neighbor,spoint,epoint))
        end
    end
    result
end
function bonds(lattice::AbstractLattice,::AcrossBonds)
    result=Bond{lattice|>keytype,lattice|>dimension}[]
    nnb=lattice|>nneighbor
    translations=reshape(product((-nnb:nnb for i=1:length(lattice.vectors))...)|>collect,:)
    for translation in translations
        remove=ntuple(i->-translation[i],length(translation))
        filter!(key->key!=remove,translations)
    end
    if length(translations)>0
        superrcoords=tile(lattice.rcoords,lattice.vectors,Tuple(translations))
        supericoords=tile(lattice.icoords,lattice.vectors,Tuple(translations))
        for link in interlinks(lattice.rcoords,superrcoords,lattice.neighbors)
            @views spoint=Point(lattice.pids[link.sindex],lattice.rcoords[:,link.sindex],lattice.icoords[:,link.sindex])
            @views epoint=Point(lattice.pids[(link.eindex-1)%size(lattice.rcoords,2)+1],superrcoords[:,link.eindex],supericoords[:,link.eindex])
            push!(result,Bond(link.neighbor,spoint,epoint))
        end
    end
    result
end

"""
    Lattice(    name::String,
                pids::Vector{<:PID},
                rcoords::AbstractMatrix{<:Real};
                icoords::AbstractMatrix{<:Real}=SMatrix{0,0,Float}(),
                vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{size(rcoords,1),Float}}(),
                neighbors::Union{Dict{Int,<:Real},Int}=1,
                coordination::Int=8
            )
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
struct Lattice{P<:PID,N} <: AbstractLattice{P,N}
    name::String
    pids::Vector{P}
    rcoords::Matrix{Float}
    icoords::Matrix{Float}
    vectors::Vector{SVector{N,Float}}
    reciprocals::Vector{SVector{N,Float}}
    neighbors::Dict{Int,Float}
    function Lattice(   name::String,
                        pids::Vector{<:PID},
                        rcoords::AbstractMatrix{<:Real};
                        icoords::AbstractMatrix{<:Real}=SMatrix{0,0,Float}(),
                        vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{size(rcoords,1),Float}}(),
                        neighbors::Union{Dict{Int,<:Real},Int}=1,
                        coordination::Int=8
                )
        length(icoords)==0 && (icoords=zeros(Float,rcoords|>size))
        @assert length(pids)==size(rcoords,2) && size(rcoords)==size(icoords) "Lattice error: dismatched shape of input pids, rcoords and icoords."
        isa(neighbors,Int) && (neighbors=Dict(i=>minlen for (i,minlen) in enumerate(minimumlengths(rcoords,vectors,neighbors,coordination=coordination))))
        new{pids|>eltype,size(rcoords,1)}(name,pids,rcoords,icoords,vectors,reciprocals(vectors),neighbors)
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
    Lattice(name,pids,rcoords,icoords=icoords,vectors=vectors,neighbors=neighbors,coordination=coordination)
end
function Lattice(   name::String,
                    sublattices::AbstractVector{<:AbstractLattice};
                    vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                    neighbors::Union{Dict{Int,<:Real},Int}=1,
                    coordination::Int=8
                    )
    @assert foldl(===,dimension.(sublattices)) "Lattice error: all sublattices should have the same dimension."
    @assert foldl(===,keytype.(sublattices)) "Lattice error: all sublattices should have the same keytype."
    len=sum(length.(sublattices))
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
    Lattice(name,pids,rcoords,icoords=icoords,vectors=vectors,neighbors=neighbors,coordination=coordination)
end

"""
    SuperLattice(   name::String,
                    sublattices::AbstractVector{<:AbstractLattice};
                    vectors::AbstractVector{<:AbstractVector{<:Real}}=SVector{0,SVector{sublattices|>eltype|>dimension,Float}}(),
                    neighbors::Dict{Int,<:Real}=Dict{Int,Float}()
                    )

SuperLattice that is composed of serveral sublattices.
"""
struct SuperLattice{L<:AbstractLattice,P<:PID,N} <: AbstractLattice{P,N}
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
        @assert foldl(===,dimension.(sublattices)) "SuperLattice error: all sublattices should have the same dimension."
        @assert foldl(===,keytype.(sublattices)) "SuperLattice error: all sublattices should have the same keytype."
        @assert all(length(sublattice.vectors)==0 for sublattice in sublattices) "SuperLattice error: all sublattices should assume open boundaries."
        len=sum(length.(sublattices))
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
        recips=reciprocals(vectors)
        new{sublattices|>eltype,sublattices|>eltype|>keytype,sublattices|>eltype|>dimension}(sublattices,name,pids,rcoords,icoords,vectors,recips,neighbors)
    end
end

struct IntraBonds <: LatticeBonds end
struct InterBonds <: LatticeBonds end
"""
    intrabonds

Indicate that bonds intra the sublattices are inquired.
!!! notes
    These bonds do not contain those accorss the periodic boundaries.
"""
const intrabonds=IntraBonds()
"""
    interbonds

Indicate that bonds inter the sublattices are inquired.
!!! notes
    These bonds do not contain those accorss the periodic boundaries.
"""
const interbonds=InterBonds()

"""
    bonds(lattice::SuperLattice,::IntraBonds) -> Vector{Bond}
    bonds(lattice::SuperLattice,::InterBonds) -> Vector{Bond}

Get the bonds of a superlattice.
"""
bonds(lattice::SuperLattice,::IntraBonds)=vcat((sublattice|>bonds for sublattice in lattice.sublattices)...)
function bonds(lattice::SuperLattice,::InterBonds)
    result=Bond{lattice|>keytype,lattice|>dimension}[]
    for (i,j) in Combinations{2}(1:length(lattice.sublattices))
        sub1,sub2=lattice.sublattices[i],lattice.sublattices[j]
        for link in interlinks(sub1.rcoords,sub2.rcoords,lattice.neighbors)
            @views spoint=Point(sub1.pids[link.sindex],sub1.rcoords[:,link.sindex],sub1.icoords[:,link.sindex])
            @views epoint=Point(sub2.pids[link.eindex],sub2.rcoords[:,link.eindex],sub2.icoords[:,link.eindex])
            push!(result,Bond(link.neighbor,spoint,epoint))
        end
    end
    result
end

"""
    Cylinder{P}(    name::String,
                    block::AbstractMatrix{<:Real},
                    translation::AbstractVector{<:Real};
                    vector::Union{AbstractVector{<:Real},Nothing}=nothing,
                    neighbors::Union{Dict{Int,<:Real},Int}=1,
                    ) where P<:PID

Cylinder of 1d and quasi 2d lattices.
"""
mutable struct Cylinder{P<:PID,N} <: AbstractLattice{P,N}
    block::Matrix{Float}
    translation::Vector{Float}
    name::String
    pids::Vector{P}
    rcoords::Matrix{Float}
    icoords::Matrix{Float}
    vectors::Vector{SVector{N,Float}}
    reciprocals::Vector{SVector{N,Float}}
    neighbors::Dict{Int,Float}
    function Cylinder{P}(   name::String,
                            block::AbstractMatrix{<:Real},
                            translation::AbstractVector{<:Real};
                            vector::Union{AbstractVector{<:Real},Nothing}=nothing,
                            neighbors::Union{Dict{Int,<:Real},Int}=1
                            ) where P<:PID
        pids=Vector{P}[]
        rcoords=zeros(Float,length(translation),0)
        icoords=zeros(Float,length(translation),0)
        vectors=vector===nothing ? SVector{length(translation),Float}[] : [vector]
        recips=reciprocals(vectors)
        isa(neighbors,Int) && (neighbors=Dict{Int,Float}(i=>Inf for i=1:neighbors))
        new{P,length(translation)}(block,translation,name,pids,rcoords,icoords,vectors,recips,neighbors)
    end
end

"""
    insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=9) where S -> Cylinder

Insert a couple of blocks into a cylinder.

The position of the cut of the cylinder is specified by the keyword argument `cut`, which is the center of the cylinder by default. All pids corresponding to a same newly inserted block share the same scope, which is specified by the parameter `ps`. Optionally, the scopes of the old pids in the cylinder can be replaced if the parameter `scopes` is assigned other than `nothing`. Note the length of `ps` is equal to the number of newly inserted blocks, while that of `scopes` should be equal to the old length of the cylinder.
"""
function Base.insert!(cylinder::Cylinder,ps::S...;cut::Int=length(cylinder)÷2+1,scopes::Union{<:AbstractVector{S},Nothing}=nothing,coordination::Int=8) where S
    @assert S===fieldtype(cylinder|>keytype,:scope) "insert! error: dismatched type of input scopes and old scopes."
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
    if any(values(cylinder.neighbors).==Inf)
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
    if any(values(neighbors).==Inf)
        minlens=minimumlengths(rcoords,vectors,cylinder|>nneighbor,coordination=coordination)
        neighbors=Dict{Int,Float}(i=>len for (i,len) in enumerate(minlens))
    end
    Lattice(name,pids,rcoords,icoords=icoords,vectors=vectors,neighbors=neighbors)
end

end #module