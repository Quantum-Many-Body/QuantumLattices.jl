module Trees

using ..Factories: Inference,Argument,Parameter,FunctionFactory,TypeFactory
using ..Factories: addfields!,addparams!,addargs!,addwhereparams!,extendbody!
using ..Factories: MixEscaped,Escaped,UnEscaped
using ..TypeTraits: efficientoperations

export treedepth,treewidth
export AbstractTree
export root,parent,children
export addnode!,deletenode!
export isleaf,level
export ancestor,descendants,siblings,leaves
export subtree,move!
export TreeCore,@tree,SimpleTree

abstract type TreeIteration end
struct TreeDepth <: TreeIteration end
struct TreeWidth <: TreeIteration end

"""
    treedepth

Indicate that the iteration over a tree is depth-first.
"""
const treedepth=TreeDepth()

"""
    treewidth

Indicate that the iteration over a tree is width-first.
"""
const treewidth=TreeWidth()

"""
    AbstractTree{Node,Data}

Abstract type for all concrete trees.
"""
abstract type AbstractTree{N,D} end

"""
    eltype(tree::AbstractTree)
    eltype(::Type{<:AbstractTree{N,D}}) where {N,D}

Get the eltype of a tree.
"""
Base.eltype(tree::AbstractTree)=tree|>typeof|>eltype
Base.eltype(::Type{<:AbstractTree{N,D}}) where {N,D}=Pair{N,D}

"""
    root(tree::AbstractTree) -> Union{keytype(tree),Nothing}

Get a tree's root node.
"""
root(tree::AbstractTree)=tree.TREECORE.root

"""
    haskey(tree::AbstractTree{N,D},node::N) where {N,D} -> Bool

Check whether a node is in a tree.
"""
Base.haskey(tree::AbstractTree{N,D},node::N) where {N,D}=haskey(tree.TREECORE.contents,node)

"""
    length(tree::AbstractTree) -> Int

Get the number of a tree's nodes.
"""
Base.length(tree::AbstractTree)=length(tree.TREECORE.contents)

"""
    parent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D} -> Union{N,Nothing}

Get the parent of a tree's node. When `node` is the tree's root, return `superparent`.
"""
parent(tree::AbstractTree{N,D},node::N,superparent::Union{N,Nothing}=nothing) where {N,D}=(node==tree|>root ? superparent : tree.TREECORE.parent[node])

"""
    children(tree::AbstractTree) -> Vector{keytype(tree)}
    children(tree::AbstractTree,::Nothing) -> Vector{keytype(tree)}
    children(tree::AbstractTree{N,D},node::N) where {N,D} -> Vector{N}

Get the children of a tree's node.
"""
children(tree::AbstractTree)=children(tree,nothing)
children(tree::AbstractTree,::Nothing)=(tree|>root===nothing ? error("children error: empty tree!") : [tree|>root])
children(tree::AbstractTree{N,D},node::N) where {N,D}=tree.TREECORE.children[node]

"""
    addnode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)
    addnode!(tree::AbstractTree{N,D},::Nothing,node::N) where {N,D} -> typeof(tree)
    addnode!(tree::AbstractTree{N,D},parent::N,node::N) where {N,D} -> typeof(tree)

Update the structure of a tree by adding a node. When the parent is `nothing`, the input tree must be empty and the input node becomes the tree's root.
"""
addnode!(tree::AbstractTree{N,D},node::N) where {N,D}=addnode!(tree,nothing,node)
function addnode!(tree::AbstractTree{N,D},::Nothing,node::N) where {N,D}
    @assert tree|>root===nothing "addnode! error: not empty tree."
    tree.TREECORE.root=node
    tree.TREECORE.children[node]=N[]
    return tree
end
function addnode!(tree::AbstractTree{N,D},parent::N,node::N) where {N,D}
    @assert haskey(tree,parent) "addnode! error: parent($parent) not in tree."
    @assert !haskey(tree,node) "addnode! error: node($node) already in tree."
    push!(tree.TREECORE.children[parent],node)
    tree.TREECORE.parent[node]=parent
    tree.TREECORE.children[node]=N[]
    return tree
end

"""
    deletenode!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)

Update the structure of a tree by deleting a node.
"""
function deletenode!(tree::AbstractTree{N,D},node::N) where {N,D}
    @assert haskey(tree,node) "deletenode! error: node($node) not in tree."
    if node==tree|>root
        tree.TREECORE.root=nothing
    else
        pnode=parent(tree,node)
        haskey(tree.TREECORE.children,pnode) && filter!(x->x!=node,tree.TREECORE.children[pnode])
    end
    pop!(tree.TREECORE.contents,node)
    pop!(tree.TREECORE.parent,node,nothing)
    pop!(tree.TREECORE.children,node)
    return tree
end

"""
    getindex(tree::AbstractTree{N,D},node::N) where {N,D} -> N

Get the data of a tree's node.
"""
Base.getindex(tree::AbstractTree{N,D},node::N) where {N,D}=tree.TREECORE.contents[node]

"""
    setindex!(tree::AbstractTree{N,D},data::D,node::N) where {N,D}

Set the data of a tree's node.
"""
Base.setindex!(tree::AbstractTree{N,D},data::D,node::N) where {N,D}=(tree.TREECORE.contents[node]=data)

"""
    empty(tree::AbstractTree)

Construct an empty tree of the same type with the input one.
"""
@generated function Base.empty(tree::AbstractTree{N,D}) where {N,D}
    contents=Expr[:(getfield(tree,$i)) for i=1:length(tree|>fieldnames)-1]
    return :(($tree)($(contents...),TreeCore{$N,$D}()))
end

"""
    keytype(tree::AbstractTree)
    keytype(::Type{<:AbstractTree{N,D}}) where {N,D}

Get a tree's node type.
"""
Base.keytype(tree::AbstractTree)=tree|>typeof|>keytype
Base.keytype(::Type{<:AbstractTree{N,D}}) where {N,D}=N

"""
    valtype(tree::AbstractTree)
    valtype(::Type{<:AbstractTree{N,D}}) where {N,D}

Get a tree's data type.
"""
Base.valtype(tree::AbstractTree)=tree|>typeof|>valtype
Base.valtype(::Type{<:AbstractTree{N,D}}) where {N,D}=D

"""
    ==(t1::T,t2::T) where T<:AbstractTree -> Bool
    isequal(t1::T,t2::T) where T<:AbstractTree -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(t1::T,t2::T) where T<:AbstractTree = ==(efficientoperations,t1,t2)
Base.isequal(t1::T,t2::T) where T<:AbstractTree=isequal(efficientoperations,t1,t2)

"""
    keys(tree::AbstractTree{N,D},::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}
    keys(tree::AbstractTree{N,D},::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}

Iterate over a tree's nodes starting from a certain `node` by depth first search or width first search.
"""
Base.keys(tree::AbstractTree{N,D},ti::TreeIteration,node::Union{N,Nothing}=tree|>root) where {N,D}=TreeKeys{typeof(ti),N,D,typeof(tree)}(tree,node)

struct TreeKeys{P<:TreeIteration,N,D,T<:AbstractTree{N,D}}
    tree::T
    node::Union{N,Nothing}
end
Base.eltype(::Type{TreeKeys{P,N,D,T}}) where {P,N,D,T}=N
Base.IteratorSize(::Type{<:TreeKeys})=Base.SizeUnknown()
Base.iterate(tk::TreeKeys)=tk.tree|>root==nothing ? nothing : (tk.node,copy(children(tk.tree,tk.node)))
function Base.iterate(tk::TreeKeys{TreeDepth,N,D,T},state::Vector{N}) where {N,D,T}
    length(state)==0 ? nothing : (node=popfirst!(state);prepend!(state,children(tk.tree,node));(node,state))
end
function Base.iterate(tk::TreeKeys{TreeWidth,N,D,T},state::Vector{N}) where {N,D,T}
    length(state)==0 ? nothing : (node=popfirst!(state);append!(state,children(tk.tree,node));(node,state))
end

"""
    values(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}
    values(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}

Iterate over a tree's data starting from a certain `node` by depth first search or width first search.
"""
Base.values(tree::AbstractTree,ti::TreeIteration,node::Union{N,Nothing}=tree|>root) where {N,D}=(tree[key] for key in keys(tree,ti,node))

"""
    pairs(tree::AbstractTree,::TreeDepth,node::Union{N,Nothing}=tree|>root) where {N,D}
    pairs(tree::AbstractTree,::TreeWidth,node::Union{N,Nothing}=tree|>root) where {N,D}

Iterate over a tree's (node,data) pairs starting from a certain `node` by depth first search or width first search.
"""
Base.pairs(tree::AbstractTree,ti::TreeIteration,node::Union{N,Nothing}=tree|>root) where {N,D}=((key,tree[key]) for key in keys(tree,ti,node))

"""
    isleaf(tree::AbstractTree{N,D},node::N) where{N,D} -> Bool

Judge whether a tree's node is a leaf (a node without children) or not.
"""
isleaf(tree::AbstractTree{N,D},node::N) where{N,D}=length(children(tree,node))==0

"""
    level(tree::AbstractTree{N,D},node::N) where {N,D} -> Int

Get the level of tree's node.
"""
function level(tree::AbstractTree{N,D},node::N) where {N,D}
    result=1
    while node!=tree|>root
        result+=1
        node=parent(tree,node)
    end
    return result
end

"""
    ancestor(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> N

Get the ancestor of a tree's node of the n-th generation.
"""
function ancestor(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D}
    @assert generation>=0 "ancestor error: generation($generation) must be non-negative."
    result=node
    for i=1:generation
        result=parent(tree,result)
    end
    result
end

"""
    descendants(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D} -> Vector{N}

Get the descendants of a tree's node of the n-th generation.
"""
function descendants(tree::AbstractTree{N,D},node::N,generation::Int=1) where {N,D}
    @assert generation>=0 "descendants error: generation($generation) must be non-negative."
    result=N[node]
    for i=1:generation
        result=vcat((children(tree,node) for node in result)...)
    end
    result
end

"""
    siblings(tree::AbstractTree{N,D},node::N) where{N,D} -> Vector{N}

Get the siblings (other nodes sharing the same parent) of a tree's node.
"""
siblings(tree::AbstractTree{N,D},node::N) where{N,D}=filter(x->x!=node,children(tree,parent(tree,node)))

"""
    leaves(tree::AbstractTree) -> Vector{keytype(tree)}

Get a tree's leaves.
"""
leaves(tree::AbstractTree)=keytype(tree)[node for node in keys(tree,treedepth) if isleaf(tree,node)]

"""
    push!(tree::AbstractTree{N,D},node::N,data::D) where {N,D} -> typeof(tree)
    push!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D} -> typeof(tree)

Push a new node to a tree. When `parent` is `nothing`, this function set the root node of an empty tree.
"""
Base.push!(tree::AbstractTree{N,D},node::N,data::D) where {N,D}=push!(tree,nothing,node,data)
function Base.push!(tree::AbstractTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D}
    addnode!(tree,parent,node)
    tree[node]=data
    return tree
end

"""
    append!(tree::AbstractTree{N,D},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)
    append!(tree::AbstractTree{N,D},node::Union{N,Nothing},subtree::AbstractTree{N,D}) where {N,D} -> typeof(tree)

Append a subtree to a tree.
"""
Base.append!(tree::AbstractTree{N,D},subtree::AbstractTree{N,D}) where {N,D}=append!(tree,nothing,subtree)
function Base.append!(tree::AbstractTree{N,D},node::Union{N,Nothing},subtree::AbstractTree{N,D}) where {N,D}
    for (key,value) in pairs(subtree,treewidth)
        push!(tree,parent(subtree,key,node),key,value)
    end
    return tree
end

"""
    delete!(tree::AbstractTree{N,D},node::N) where {N,D} -> typeof(tree)

Delete a node and all its descendants from a tree.
"""
function Base.delete!(tree::AbstractTree{N,D},node::N) where {N,D}
    for key in collect(N,keys(tree,treedepth,node))
        deletenode!(tree,key)
    end
    return tree
end

"""
    empty!(tree::AbstractTree) -> typeof(tree)

Empty a tree.
"""
Base.empty!(tree::AbstractTree)=delete!(tree,tree|>root)

"""
    subtree(tree::AbstractTree{N,D},node::N) where{N,D} -> typeof(tree)

Get a subtree whose root is `node`.
"""
function subtree(tree::AbstractTree{N,D},node::N) where{N,D}
    result=empty(tree)
    for (i,(key,value)) in enumerate(pairs(tree,treedepth,node))
        push!(result,i==1 ? nothing : parent(tree,key),key,value)
    end
    return result
end

"""
    move!(tree::AbstractTree{N,D},node::N,parent::N) where {N,D} -> typeof(tree)

Move a subtree to a new position.
"""
function move!(tree::AbstractTree{N,D},node::N,parent::N) where {N,D}
    sub=subtree(tree,node)
    delete!(tree,node)
    append!(tree,parent,sub)
    return tree
end

"""
    TreeCore()

The core of a tree.
"""
mutable struct TreeCore{N,D}
    root::Union{N,Nothing}
    contents::Dict{N,D}
    parent::Dict{N,N}
    children::Dict{N,Vector{N}}
    TreeCore{N,D}() where {N,D}=new{N,D}(nothing,Dict{N,D}(),Dict{N,N}(),Dict{N,Vector{N}}())
end

"""
    ==(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool
    isequal(tc1::TC,tc2::TC) where TC<:TreeCore -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(tc1::TC,tc2::TC) where TC<:TreeCore = ==(efficientoperations,tc1,tc2)
Base.isequal(tc1::TC,tc2::TC) where TC<:TreeCore=isequal(efficientoperations,tc1,tc2)

"""
    @tree structdef treeparams::Union{Expr,Nothing}=nothing

Decorate a "raw" struct to be a subtype of `AbstractTree`.
!!! note
    1. A "raw" struct means:
       - It has no explicit supertype;
       - It has no inner constructor;
       - It has no attribute `:TREECORE`.
    2. The keytype and valtype can be assigned by the argument `treeparams` in the form `{keytype,valtype}`.
       - When the formal argument names of keytype and valtype are not assigned, they can be automatically generated by the functioin `gensym`.
         For example, all of the structs after the decration by the following codes
         ```julia
         @tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end)
         @tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,
               {::String,::Int}
               )
         @tree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,
               {<:AbstractString,<:Number}
               )
         ```
         will have three type parameters.
       - When the formal argument names of keytype and valtype overlap with those of the raw struct type parameters, the duplicates will be considered as the same.
         For example, the decorated struct `SubTreeWithOverlappedParametricFields` by the following code
         ```julia
         @tree(struct TreeWithOverlappedParametricFields{N} info::Vector{N} end,
               {N<:AbstractString,D<:Number}
               )
         ```
         only has two type parameters `N<:AbstractString` and `D<:Number`, where the `N` in the `info::Vector{N}` is the same `N` with that in the decorated attribute `TREECORE::TreeCore{N,D}`.
       - When the formal argument names of keytype and valtype have no intersection with those of the raw struct type parameters,
         the type parameters of the decorated struct will be just extended by keytype and valtype.
         For example, the decorated struct `SubTreeWithParametricFields` by the following code
         ```julia
         @tree(struct TreeWithParametricFields{T} info::Vector{T} end,
               {N<:AbstractString,D<:Number}
               )
         ```
         have 3 type parameters, `T`, `N<:AbstractString` and `D<:Number`.
"""
macro tree(structdef,treeparams::Union{Expr,Nothing}=nothing)
    tf=TypeFactory(structdef)
    fieldnames=[field.name for field in tf.fields]
    paramnames=[param.name for param in tf.params]
    @assert tf.supertype==Inference(:Any) "@tree error: no explicit supertype except `Any` is allowed."
    @assert length(tf.constructors)==0 "@tree error: no inner constructor is allowed."
    @assert :TREECORE ∉ fieldnames "@tree error: :TREECORE is a reserved attribute name."
    if treeparams==nothing
        treeparamnames=(gensym(),gensym())
        treeparams=[Parameter(treeparamnames[1]),Parameter(treeparamnames[2])]
    else
        @assert (treeparams.head==:braces && length(treeparams.args)==2) "@tree error: not supported treeparams."
        treeparams=[Parameter(arg) for arg in treeparams.args]
        treeparams[1].name===nothing && (treeparams[1].name=gensym())
        treeparams[2].name===nothing && (treeparams[2].name=gensym())
        treeparamnames=tuple((tp.name for tp in treeparams)...)
        filter!(tp->tp.name ∉ paramnames,treeparams)
    end
    tf.supertype=Inference(:(AbstractTree{$(treeparamnames...)}))
    append!(paramnames,(tp.name for tp in treeparams))
    addfields!(tf,:(TREECORE::TreeCore{$(treeparamnames...)}))
    addparams!(tf,treeparams...)
    outer=FunctionFactory(name=tf.name)
    addparams!(outer,treeparamnames...)
    addargs!(outer,(Argument(name=field.name,type=field.type) for field in tf.fields[1:end-1])...)
    addwhereparams!(outer,tf.params...)
    extendbody!(outer,:($(tf.name)($(fieldnames...),TreeCore{$(treeparamnames...)}())))
    sbtreedef=tf(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...,:AbstractTree,:TreeCore)))
    outerdef=outer(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...)))
    return Expr(:block,:(Base.@__doc__($sbtreedef)),outerdef)
end

"""
    SimpleTree{N,D}() where {N,D}

The minimum tree structure that implements all the default tree methods.
"""
struct SimpleTree{N,D} <: AbstractTree{N,D}
    TREECORE::TreeCore{N,D}
end
SimpleTree{N,D}() where {N,D}=SimpleTree(TreeCore{N,D}())

end #module
