module Tree

const treedepth=Val(0)
const treewidth=Val(1)

"""
    AbstractTree{Node,Data}

Abstract type for all concrete trees.
"""
abstract type AbstractTree{N,D} end

Base.eltype(::Type{<:AbstractTree{N,D}}) where {N,D}=Pair{N,D}
Base.eltype(at::AbstractTree)=at|>typeof|>eltype

Base.keytype(::Type{<:AbstractTree{N,D}}) where {N,D}=N
Base.keytype(at::AbstractTree)=at|>typeof|>keytype

Base.valtype(::Type{<:AbstractTree{N,D}}) where {N,D}=D
Base.valtype(at::AbstractTree)=at|>typeof|>valtype

Base.push!(at::AbstractTree{N,D},parent::N,node::N,data::D) where {N,D}=nothing

Base.append!(at::AbstractTree{N,D},parent::N,subtree::AbstractTree{N,D}) where {N,D}=nothing

Base.delete!(at::AbstractTree{N,D},parent::N) where {N,D}=nothing

Base.keys(at::AbstractTree,::typeof(treedepth))=nothing
Base.keys(at::AbstractTree,::typeof(treewidth))=nothing

Base.values(at::AbstractTree,::typeof(treedepth))=nothing
Base.values(at::AbstractTree,::typeof(treewidth))=nothing

Base.pairs(at::AbstractTree,::typeof(treedepth))=nothing
Base.pairs(at::AbstractTree,::typeof(treewidth))=nothing

Base.empty!(at::AbstractTree)=nothing

move(at::AbstractTree{N,D},node::N,parent::N) where {N,D}=nothing

parent(at::AbstractTree{N,D},node::N) where {N,D}=nothing

children(at::AbstractTree{N,D},node::N) where {N,D}=nothing

ancestor(at::AbstractTree{N,D},node::N,generation::Int=1) where {N,D}=nothing

descendants(at::AbstractTree{N,D},node::N,generation::Int=1) where {N,D}=nothing

siblings(at::AbstractTree{N,D},node::N) where{N,D}=nothing

subtree(at::AbstractTree{N,D},node::N) where{N,D}=nothing

isleaf(at::AbstractTree{N,D},node::N) where{N,D}=nothing

leaves(at::AbstractTree)=nothing

level(at::AbstractTree{N,D},node::N) where {N,D}=nothing

end #module
