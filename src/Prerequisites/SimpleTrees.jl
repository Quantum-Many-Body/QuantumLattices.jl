module SimpleTrees

using ..Factories: Inference,Argument,Parameter,FunctionFactory,TypeFactory
using ..Factories: addfields!,addparams!,addargs!,addwhereparams!,extendbody!
using ..Factories: MixEscaped,Escaped,UnEscaped
using ..TypeTraits: efficientoperations

export simpletreedepth,simpletreewidth
export AbstractSimpleTree
export root,parent,children
export addnode!,deletenode!
export isleaf,level
export ancestor,descendants,siblings,leaves
export subtree,move!
export SimpleTreeCore,@simpletree,SimpleTree

abstract type SimpleTreeIteration end
struct SimpleTreeDepth <: SimpleTreeIteration end
struct SimpleTreeWidth <: SimpleTreeIteration end

"""
    simpletreedepth

Indicate that the iteration over a tree is depth-first.
"""
const simpletreedepth=SimpleTreeDepth()

"""
    simpletreewidth

Indicate that the iteration over a tree is width-first.
"""
const simpletreewidth=SimpleTreeWidth()

"""
    AbstractSimpleTree{Node,Data}

Abstract type for all concrete trees.
"""
abstract type AbstractSimpleTree{N,D} end

"""
    eltype(tree::AbstractSimpleTree)
    eltype(::Type{<:AbstractSimpleTree{N,D}}) where {N,D}

Get the eltype of a tree.
"""
Base.eltype(tree::AbstractSimpleTree)=tree|>typeof|>eltype
Base.eltype(::Type{<:AbstractSimpleTree{N,D}}) where {N,D}=Pair{N,D}

"""
    root(tree::AbstractSimpleTree) -> Union{keytype(tree),Nothing}

Get a tree's root node.
"""
root(tree::AbstractSimpleTree)=tree.TREECORE.root

"""
    haskey(tree::AbstractSimpleTree{N},node::N) where N -> Bool

Check whether a node is in a tree.
"""
Base.haskey(tree::AbstractSimpleTree{N},node::N) where N=haskey(tree.TREECORE.contents,node)

"""
    length(tree::AbstractSimpleTree) -> Int

Get the number of a tree's nodes.
"""
Base.length(tree::AbstractSimpleTree)=length(tree.TREECORE.contents)

"""
    parent(tree::AbstractSimpleTree{N},node::N,superparent::Union{N,Nothing}=nothing) where N -> Union{N,Nothing}

Get the parent of a tree's node. When `node` is the tree's root, return `superparent`.
"""
parent(tree::AbstractSimpleTree{N},node::N,superparent::Union{N,Nothing}=nothing) where N=(node==tree|>root ? superparent : tree.TREECORE.parent[node])

"""
    children(tree::AbstractSimpleTree) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree,::Nothing) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree{N},node::N) where N -> Vector{N}

Get the children of a tree's node.
"""
children(tree::AbstractSimpleTree)=children(tree,nothing)
children(tree::AbstractSimpleTree,::Nothing)=(tree|>root===nothing ? error("children error: empty tree!") : [tree|>root])
children(tree::AbstractSimpleTree{N},node::N) where N=tree.TREECORE.children[node]

"""
    addnode!(tree::AbstractSimpleTree{N},node::N) where N} -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N},::Nothing,node::N) where N -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N},parent::N,node::N) where N -> typeof(tree)

Update the structure of a tree by adding a node. When the parent is `nothing`, the input tree must be empty and the input node becomes the tree's root.
"""
addnode!(tree::AbstractSimpleTree{N},node::N) where N=addnode!(tree,nothing,node)
function addnode!(tree::AbstractSimpleTree{N},::Nothing,node::N) where N
    @assert tree|>root===nothing "addnode! error: not empty tree."
    tree.TREECORE.root=node
    tree.TREECORE.children[node]=N[]
    return tree
end
function addnode!(tree::AbstractSimpleTree{N},parent::N,node::N) where N
    @assert haskey(tree,parent) "addnode! error: parent($parent) not in tree."
    @assert !haskey(tree,node) "addnode! error: node($node) already in tree."
    push!(tree.TREECORE.children[parent],node)
    tree.TREECORE.parent[node]=parent
    tree.TREECORE.children[node]=N[]
    return tree
end

"""
    deletenode!(tree::AbstractSimpleTree{N},node::N) where N -> typeof(tree)

Update the structure of a tree by deleting a node.
"""
function deletenode!(tree::AbstractSimpleTree{N},node::N) where N
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
    getindex(tree::AbstractSimpleTree{N},node::N) where N -> N

Get the data of a tree's node.
"""
Base.getindex(tree::AbstractSimpleTree{N},node::N) where N=tree.TREECORE.contents[node]

"""
    setindex!(tree::AbstractSimpleTree{N,D},data::D,node::N) where {N,D}

Set the data of a tree's node.
"""
Base.setindex!(tree::AbstractSimpleTree{N,D},data::D,node::N) where {N,D}=(tree.TREECORE.contents[node]=data)

"""
    empty(tree::AbstractSimpleTree)

Construct an empty tree of the same type with the input one.
"""
@generated function Base.empty(tree::AbstractSimpleTree{N,D}) where {N,D}
    contents=Expr[:(getfield(tree,$i)) for i=1:length(tree|>fieldnames)-1]
    return :(($tree)($(contents...),SimpleTreeCore{$N,$D}()))
end

"""
    keytype(tree::AbstractSimpleTree)
    keytype(::Type{<:AbstractSimpleTree{N}}) where N

Get a tree's node type.
"""
Base.keytype(tree::AbstractSimpleTree)=tree|>typeof|>keytype
Base.keytype(::Type{<:AbstractSimpleTree{N}}) where N=N

"""
    valtype(tree::AbstractSimpleTree)
    valtype(::Type{<:AbstractSimpleTree{N,D} where N}) where D=D

Get a tree's data type.
"""
Base.valtype(tree::AbstractSimpleTree)=tree|>typeof|>valtype
Base.valtype(::Type{<:AbstractSimpleTree{N,D} where N}) where D=D

"""
    ==(t1::T,t2::T) where T<:AbstractSimpleTree -> Bool
    isequal(t1::T,t2::T) where T<:AbstractSimpleTree -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(t1::T,t2::T) where T<:AbstractSimpleTree = ==(efficientoperations,t1,t2)
Base.isequal(t1::T,t2::T) where T<:AbstractSimpleTree=isequal(efficientoperations,t1,t2)

"""
    keys(tree::AbstractSimpleTree{N},::SimpleTreeDepth,node::Union{N,Nothing}=tree|>root) where N
    keys(tree::AbstractSimpleTree{N},::SimpleTreeWidth,node::Union{N,Nothing}=tree|>root) where N

Iterate over a tree's nodes starting from a certain `node` by depth first search or width first search.
"""
Base.keys(tree::AbstractSimpleTree{N},ti::SimpleTreeIteration,node::Union{N,Nothing}=tree|>root) where N=TreeKeys{typeof(ti),N,valtype(tree),typeof(tree)}(tree,node)
struct TreeKeys{P<:SimpleTreeIteration,N,D,T<:AbstractSimpleTree{N,D}}
    tree::T
    node::Union{N,Nothing}
end
Base.eltype(::Type{TreeKeys{<:SimpleTreeIteration,N,D,<:AbstractSimpleTree{N,D}} where D}) where N=N
Base.IteratorSize(::Type{<:TreeKeys})=Base.SizeUnknown()
Base.iterate(tk::TreeKeys)=tk.tree|>root==nothing ? nothing : (tk.node,copy(children(tk.tree,tk.node)))
function Base.iterate(tk::TreeKeys{SimpleTreeDepth,N},state::Vector{N}) where N
    length(state)==0 ? nothing : (node=popfirst!(state);prepend!(state,children(tk.tree,node));(node,state))
end
function Base.iterate(tk::TreeKeys{SimpleTreeWidth,N},state::Vector{N}) where N
    length(state)==0 ? nothing : (node=popfirst!(state);append!(state,children(tk.tree,node));(node,state))
end

"""
    values(tree::AbstractSimpleTree{N},::SimpleTreeDepth,node::Union{N,Nothing}=tree|>root) where N
    values(tree::AbstractSimpleTree{N},::SimpleTreeWidth,node::Union{N,Nothing}=tree|>root) where N

Iterate over a tree's data starting from a certain `node` by depth first search or width first search.
"""
Base.values(tree::AbstractSimpleTree{N},ti::SimpleTreeIteration,node::Union{N,Nothing}=tree|>root) where N=(tree[key] for key in keys(tree,ti,node))

"""
    pairs(tree::AbstractSimpleTree{N},::SimpleTreeDepth,node::Union{N,Nothing}=tree|>root) where N
    pairs(tree::AbstractSimpleTree{N},::SimpleTreeWidth,node::Union{N,Nothing}=tree|>root) where N

Iterate over a tree's (node,data) pairs starting from a certain `node` by depth first search or width first search.
"""
Base.pairs(tree::AbstractSimpleTree{N},ti::SimpleTreeIteration,node::Union{N,Nothing}=tree|>root) where N=((key,tree[key]) for key in keys(tree,ti,node))

"""
    isleaf(tree::AbstractSimpleTree{N},node::N) where N -> Bool

Judge whether a tree's node is a leaf (a node without children) or not.
"""
isleaf(tree::AbstractSimpleTree{N},node::N) where N=length(children(tree,node))==0

"""
    level(tree::AbstractSimpleTree{N},node::N) where N -> Int

Get the level of tree's node.
"""
function level(tree::AbstractSimpleTree{N},node::N) where N
    result=1
    while node!=tree|>root
        result+=1
        node=parent(tree,node)
    end
    return result
end

"""
    ancestor(tree::AbstractSimpleTree{N},node::N,generation::Int=1) where N -> N

Get the ancestor of a tree's node of the n-th generation.
"""
function ancestor(tree::AbstractSimpleTree{N},node::N,generation::Int=1) where N
    @assert generation>=0 "ancestor error: generation($generation) must be non-negative."
    result=node
    for i=1:generation
        result=parent(tree,result)
    end
    result
end

"""
    descendants(tree::AbstractSimpleTree{N},node::N,generation::Int=1) where N -> Vector{N}

Get the descendants of a tree's node of the n-th generation.
"""
function descendants(tree::AbstractSimpleTree{N},node::N,generation::Int=1) where N
    @assert generation>=0 "descendants error: generation($generation) must be non-negative."
    result=N[node]
    for i=1:generation
        result=vcat((children(tree,node) for node in result)...)
    end
    result
end

"""
    siblings(tree::AbstractSimpleTree{N},node::N) where N -> Vector{N}

Get the siblings (other nodes sharing the same parent) of a tree's node.
"""
siblings(tree::AbstractSimpleTree{N},node::N) where N=filter(x->x!=node,children(tree,parent(tree,node)))

"""
    leaves(tree::AbstractSimpleTree) -> Vector{keytype(tree)}

Get a tree's leaves.
"""
leaves(tree::AbstractSimpleTree)=keytype(tree)[node for node in keys(tree,simpletreedepth) if isleaf(tree,node)]

"""
    push!(tree::AbstractSimpleTree{N,D},node::N,data::D) where {N,D} -> typeof(tree)
    push!(tree::AbstractSimpleTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D} -> typeof(tree)

Push a new node to a tree. When `parent` is `nothing`, this function set the root node of an empty tree.
"""
Base.push!(tree::AbstractSimpleTree{N,D},node::N,data::D) where {N,D}=push!(tree,nothing,node,data)
function Base.push!(tree::AbstractSimpleTree{N,D},parent::Union{N,Nothing},node::N,data::D) where {N,D}
    addnode!(tree,parent,node)
    tree[node]=data
    return tree
end

"""
    append!(tree::AbstractSimpleTree{N,D},subtree::AbstractSimpleTree{N,D}) where {N,D} -> typeof(tree)
    append!(tree::AbstractSimpleTree{N,D},node::Union{N,Nothing},subtree::AbstractSimpleTree{N,D}) where {N,D} -> typeof(tree)

Append a subtree to a tree.
"""
Base.append!(tree::AbstractSimpleTree{N,D},subtree::AbstractSimpleTree{N,D}) where {N,D}=append!(tree,nothing,subtree)
function Base.append!(tree::AbstractSimpleTree{N,D},node::Union{N,Nothing},subtree::AbstractSimpleTree{N,D}) where {N,D}
    for (key,value) in pairs(subtree,simpletreewidth)
        push!(tree,parent(subtree,key,node),key,value)
    end
    return tree
end

"""
    delete!(tree::AbstractSimpleTree{N},node::N) where N -> typeof(tree)

Delete a node and all its descendants from a tree.
"""
function Base.delete!(tree::AbstractSimpleTree{N},node::N) where N
    for key in collect(N,keys(tree,simpletreedepth,node))
        deletenode!(tree,key)
    end
    return tree
end

"""
    empty!(tree::AbstractSimpleTree) -> typeof(tree)

Empty a tree.
"""
Base.empty!(tree::AbstractSimpleTree)=delete!(tree,tree|>root)

"""
    subtree(tree::AbstractSimpleTree{N},node::N) where N -> typeof(tree)

Get a subtree whose root is `node`.
"""
function subtree(tree::AbstractSimpleTree{N},node::N) where N
    result=empty(tree)
    for (i,(key,value)) in enumerate(pairs(tree,simpletreedepth,node))
        push!(result,i==1 ? nothing : parent(tree,key),key,value)
    end
    return result
end

"""
    move!(tree::AbstractSimpleTree{N},node::N,parent::N) where N -> typeof(tree)

Move a subtree to a new position.
"""
function move!(tree::AbstractSimpleTree{N},node::N,parent::N) where N
    sub=subtree(tree,node)
    delete!(tree,node)
    append!(tree,parent,sub)
    return tree
end

"""
    SimpleTreeCore()

The core of a tree.
"""
mutable struct SimpleTreeCore{N,D}
    root::Union{N,Nothing}
    contents::Dict{N,D}
    parent::Dict{N,N}
    children::Dict{N,Vector{N}}
    SimpleTreeCore{N,D}() where {N,D}=new{N,D}(nothing,Dict{N,D}(),Dict{N,N}(),Dict{N,Vector{N}}())
end

"""
    ==(tc1::TC,tc2::TC) where TC<:SimpleTreeCore -> Bool
    isequal(tc1::TC,tc2::TC) where TC<:SimpleTreeCore -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(tc1::TC,tc2::TC) where TC<:SimpleTreeCore = ==(efficientoperations,tc1,tc2)
Base.isequal(tc1::TC,tc2::TC) where TC<:SimpleTreeCore=isequal(efficientoperations,tc1,tc2)

"""
    @simpletree structdef treeparams::Union{Expr,Nothing}=nothing

Decorate a "raw" struct to be a subtype of `AbstractSimpleTree`.
!!! note
    1. A "raw" struct means:
       - It has no explicit supertype;
       - It has no inner constructor;
       - It has no attribute `:TREECORE`.
    2. The keytype and valtype can be assigned by the argument `treeparams` in the form `{keytype,valtype}`.
       - When the formal argument names of keytype and valtype are not assigned, they can be automatically generated by the functioin `gensym`.
         For example, all of the structs after the decration by the following codes
         ```julia
         @simpletree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end)
         @simpletree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,
               {::String,::Int}
               )
         @simpletree(struct TreeWithWrongTypeParameterNames{N} info::Vector{N} end,
               {<:AbstractString,<:Number}
               )
         ```
         will have three type parameters.
       - When the formal argument names of keytype and valtype overlap with those of the raw struct type parameters, the duplicates will be considered as the same.
         For example, the decorated struct `SubTreeWithOverlappedParametricFields` by the following code
         ```julia
         @simpletree(struct TreeWithOverlappedParametricFields{N} info::Vector{N} end,
               {N<:AbstractString,D<:Number}
               )
         ```
         only has two type parameters `N<:AbstractString` and `D<:Number`, where the `N` in the `info::Vector{N}` is the same `N` with that in the decorated attribute `TREECORE::SimpleTreeCore{N,D}`.
       - When the formal argument names of keytype and valtype have no intersection with those of the raw struct type parameters,
         the type parameters of the decorated struct will be just extended by keytype and valtype.
         For example, the decorated struct `SubTreeWithParametricFields` by the following code
         ```julia
         @simpletree(struct TreeWithParametricFields{T} info::Vector{T} end,
               {N<:AbstractString,D<:Number}
               )
         ```
         have 3 type parameters, `T`, `N<:AbstractString` and `D<:Number`.
"""
macro simpletree(structdef,treeparams::Union{Expr,Nothing}=nothing)
    tf=TypeFactory(structdef)
    fieldnames=[field.name for field in tf.fields]
    paramnames=[param.name for param in tf.params]
    @assert tf.supertype==Inference(:Any) "@simpletree error: no explicit supertype except `Any` is allowed."
    @assert length(tf.constructors)==0 "@simpletree error: no inner constructor is allowed."
    @assert :TREECORE ∉ fieldnames "@simpletree error: :TREECORE is a reserved attribute name."
    if treeparams==nothing
        treeparamnames=(gensym(),gensym())
        treeparams=[Parameter(treeparamnames[1]),Parameter(treeparamnames[2])]
    else
        @assert (treeparams.head==:braces && length(treeparams.args)==2) "@simpletree error: not supported treeparams."
        treeparams=[Parameter(arg) for arg in treeparams.args]
        treeparams[1].name===nothing && (treeparams[1].name=gensym())
        treeparams[2].name===nothing && (treeparams[2].name=gensym())
        treeparamnames=tuple((tp.name for tp in treeparams)...)
        filter!(tp->tp.name ∉ paramnames,treeparams)
    end
    tf.supertype=Inference(:(AbstractSimpleTree{$(treeparamnames...)}))
    append!(paramnames,(tp.name for tp in treeparams))
    addfields!(tf,:(TREECORE::SimpleTreeCore{$(treeparamnames...)}))
    addparams!(tf,treeparams...)
    outer=FunctionFactory(name=tf.name)
    addparams!(outer,treeparamnames...)
    addargs!(outer,(Argument(name=field.name,type=field.type) for field in tf.fields[1:end-1])...)
    addwhereparams!(outer,tf.params...)
    extendbody!(outer,:($(tf.name)($(fieldnames...),SimpleTreeCore{$(treeparamnames...)}())))
    sbtreedef=tf(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...,:AbstractSimpleTree,:SimpleTreeCore)))
    outerdef=outer(MixEscaped(Escaped(tf.name),UnEscaped(paramnames...)))
    return Expr(:block,:(Base.@__doc__($sbtreedef)),outerdef)
end

"""
    SimpleTree{N,D}() where {N,D}

The minimum tree structure that implements all the default tree methods.
"""
struct SimpleTree{N,D} <: AbstractSimpleTree{N,D}
    TREECORE::SimpleTreeCore{N,D}
end
SimpleTree{N,D}() where {N,D}=SimpleTree(SimpleTreeCore{N,D}())

end #module
