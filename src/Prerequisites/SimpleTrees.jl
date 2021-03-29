module SimpleTrees

using ..Traits: getcontent, parametertype, efficientoperations, rawtype

import ..Traits: contentnames, dissolve

export simpletreedepth, simpletreewidth
export SimpleTreeCore, AbstractSimpleTree, SimpleTree
export root, parent, children, ancestor, descendants, siblings, leaves, subtree
export isleaf, level, addnode!, deletenode!, move!

"""
    SimpleTreeCore()

The core of a simple tree.
"""
mutable struct SimpleTreeCore{N, D}
    root::Union{N, Nothing}
    contents::Dict{N, D}
    parent::Dict{N, N}
    children::Dict{N, Vector{N}}
    SimpleTreeCore{N, D}() where {N, D} = new{N, D}(nothing, Dict{N, D}(), Dict{N, N}(), Dict{N, Vector{N}}())
end

"""
    ==(tc1::TC, tc2::TC) where TC<:SimpleTreeCore -> Bool
    isequal(tc1::TC, tc2::TC) where TC<:SimpleTreeCore -> Bool

Overloaded equivalent operator.
"""
@inline Base.:(==)(tc1::TC, tc2::TC) where {TC<:SimpleTreeCore} = ==(efficientoperations, tc1, tc2)
@inline Base.isequal(tc1::TC, tc2::TC) where {TC<:SimpleTreeCore} = isequal(efficientoperations, tc1, tc2)

abstract type SimpleTreeIteration end
struct SimpleTreeDepth <: SimpleTreeIteration end
struct SimpleTreeWidth <: SimpleTreeIteration end

"""
    simpletreedepth

Indicate that the iteration over a tree is depth-first.
"""
const simpletreedepth = SimpleTreeDepth()

"""
    simpletreewidth

Indicate that the iteration over a tree is width-first.
"""
const simpletreewidth = SimpleTreeWidth()

"""
    AbstractSimpleTree{Node, Data}

Abstract type for all concrete trees.
"""
abstract type AbstractSimpleTree{N, D} end
@inline contentnames(::Type{<:AbstractSimpleTree}) = (:TREECORE,)

"""
    SimpleTreeCore(tree::AbstractSimpleTree) -> SimpleTreeCore

Get the core of a simple tree.
"""
@inline SimpleTreeCore(tree::AbstractSimpleTree) = getcontent(tree, :TREECORE)

"""
    keytype(tree::AbstractSimpleTree)
    keytype(::Type{T}) where {T<:AbstractSimpleTree}

Get a tree's node type.
"""
@inline Base.keytype(tree::AbstractSimpleTree) = keytype(typeof(tree))
@inline @generated Base.keytype(::Type{T}) where {T<:AbstractSimpleTree} = parametertype(supertype(T, :AbstractSimpleTree), 1)

"""
    valtype(tree::AbstractSimpleTree)
    valtype(::Type{T}) where {T<:AbstractSimpleTree}

Get a tree's data type.
"""
@inline Base.valtype(tree::AbstractSimpleTree) = valtype(typeof(tree))
@inline @generated Base.valtype(::Type{T}) where {T<:AbstractSimpleTree} = parametertype(supertype(T, :AbstractSimpleTree), 2)

"""
    eltype(tree::AbstractSimpleTree)
    eltype(::Type{T}) where {T<:AbstractSimpleTree}

Get the eltype of a tree.
"""
@inline Base.eltype(tree::AbstractSimpleTree) = eltype(typeof(tree))
@inline Base.eltype(::Type{T}) where {T<:AbstractSimpleTree} = Pair{keytype(T), valtype(T)}

"""
    root(tree::AbstractSimpleTree) -> Union{keytype(tree), Nothing}

Get a tree's root node.
"""
@inline root(tree::AbstractSimpleTree) = SimpleTreeCore(tree).root

"""
    haskey(tree::AbstractSimpleTree{N}, node::N) where N -> Bool

Check whether a node is in a tree.
"""
@inline Base.haskey(tree::AbstractSimpleTree{N}, node::N) where N = haskey(SimpleTreeCore(tree).contents, node)

"""
    length(tree::AbstractSimpleTree) -> Int

Get the number of a tree's nodes.
"""
@inline Base.length(tree::AbstractSimpleTree) = length(SimpleTreeCore(tree).contents)

"""
    parent(tree::AbstractSimpleTree{N}, node::N, superparent::Union{N, Nothing}=nothing) where N -> Union{N, Nothing}

Get the parent of a tree's node. When `node` is the tree's root, return `superparent`.
"""
@inline parent(tree::AbstractSimpleTree{N}, node::N, superparent::Union{N, Nothing}=nothing) where {N} = (node == root(tree) ? superparent : SimpleTreeCore(tree).parent[node])

"""
    children(tree::AbstractSimpleTree) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree, ::Nothing) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree{N}, node::N) where N -> Vector{N}

Get the children of a tree's node.
"""
@inline children(tree::AbstractSimpleTree) = children(tree, nothing)
@inline children(tree::AbstractSimpleTree, ::Nothing) = (root(tree) === nothing) ? error("children error: empty tree!") : [root(tree)]
@inline children(tree::AbstractSimpleTree{N}, node::N) where N = SimpleTreeCore(tree).children[node]

"""
    addnode!(tree::AbstractSimpleTree{N}, node::N) where N} -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N}, ::Nothing, node::N) where N -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N}, parent::N, node::N) where N -> typeof(tree)

Update the structure of a tree by adding a node. When the parent is `nothing`, the input tree must be empty and the input node becomes the tree's root.
"""
@inline addnode!(tree::AbstractSimpleTree{N}, node::N) where {N} = addnode!(tree, nothing, node)
function addnode!(tree::AbstractSimpleTree{N}, ::Nothing, node::N) where N
    @assert root(tree) === nothing "addnode! error: not empty tree."
    SimpleTreeCore(tree).root = node
    SimpleTreeCore(tree).children[node] = N[]
    return tree
end
function addnode!(tree::AbstractSimpleTree{N}, parent::N, node::N) where N
    @assert haskey(tree, parent) "addnode! error: parent($parent) not in tree."
    @assert !haskey(tree, node) "addnode! error: node($node) already in tree."
    push!(SimpleTreeCore(tree).children[parent], node)
    SimpleTreeCore(tree).parent[node] = parent
    SimpleTreeCore(tree).children[node] = N[]
    return tree
end

"""
    deletenode!(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Update the structure of a tree by deleting a node.
"""
function deletenode!(tree::AbstractSimpleTree{N}, node::N) where N
    @assert haskey(tree, node) "deletenode! error: node($node) not in tree."
    if node == root(tree)
        SimpleTreeCore(tree).root = nothing
    else
        pnode = parent(tree, node)
        haskey(SimpleTreeCore(tree).children, pnode) && filter!(!=(node), SimpleTreeCore(tree).children[pnode])
    end
    pop!(SimpleTreeCore(tree).contents, node)
    pop!(SimpleTreeCore(tree).parent, node, nothing)
    pop!(SimpleTreeCore(tree).children, node)
    return tree
end

"""
    getindex(tree::AbstractSimpleTree{N}, node::N) where N -> N

Get the data of a tree's node.
"""
@inline Base.getindex(tree::AbstractSimpleTree{N}, node::N) where {N} = SimpleTreeCore(tree).contents[node]

"""
    setindex!(tree::AbstractSimpleTree{N, D}, data::D, node::N) where {N, D}

Set the data of a tree's node.
"""
@inline Base.setindex!(tree::AbstractSimpleTree{N, D}, data::D, node::N) where {N, D} = (SimpleTreeCore(tree).contents[node] = data)

"""
    empty(tree::AbstractSimpleTree)

Construct an empty tree of the same type with the input one.
"""
@inline Base.empty(tree::AbstractSimpleTree) = rawtype(typeof(tree))(dissolve(tree, empty)...)
@inline dissolve(tree::AbstractSimpleTree{N, D}, ::Val{:TREECORE}, ::typeof(empty), args::Tuple, kwargs::NamedTuple) where {N, D} = SimpleTreeCore{N, D}()

"""
    ==(t1::T, t2::T) where {T<:AbstractSimpleTree} -> Bool

Overloaded equivalent operator.
"""
@inline Base.:(==)(t1::T, t2::T) where {T<:AbstractSimpleTree} = ==(efficientoperations, t1, t2)

"""
    isequal(t1::T, t2::T) where {T<:AbstractSimpleTree} -> Bool

Overloaded equivalent operator.
"""
@inline Base.isequal(t1::T, t2::T) where {T<:AbstractSimpleTree} = isequal(efficientoperations, t1, t2)

"""
    keys(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    keys(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's nodes starting from a certain `node` by depth first search or width first search.
"""
function Base.keys(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where N
    TreeKeys{typeof(ti), N, valtype(tree), typeof(tree)}(tree, node)
end
struct TreeKeys{P<:SimpleTreeIteration, N, D, T<:AbstractSimpleTree{N, D}}
    tree::T
    node::Union{N, Nothing}
end
Base.eltype(::Type{<:TreeKeys{<:SimpleTreeIteration, N, D, <:AbstractSimpleTree{N, D}} where D}) where {N} = N
Base.IteratorSize(::Type{<:TreeKeys}) = Base.SizeUnknown()
Base.iterate(tk::TreeKeys) = (root(tk.tree) === nothing) ? nothing : (tk.node, copy(children(tk.tree, tk.node)))
function Base.iterate(tk::TreeKeys{SimpleTreeDepth, N}, state::Vector{N}) where N
    (length(state) == 0) ? nothing : (node = popfirst!(state); prepend!(state, children(tk.tree, node)); (node, state))
end
function Base.iterate(tk::TreeKeys{SimpleTreeWidth, N}, state::Vector{N}) where N
    (length(state) == 0) ? nothing : (node = popfirst!(state); append!(state, children(tk.tree, node)); (node, state))
end

"""
    values(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    values(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's data starting from a certain `node` by depth first search or width first search.
"""
@inline Base.values(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where {N} = (tree[key] for key in keys(tree, ti, node))

"""
    pairs(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    pairs(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's (node, data) pairs starting from a certain `node` by depth first search or width first search.
"""
@inline Base.pairs(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where {N} = ((key, tree[key]) for key in keys(tree, ti, node))

"""
    isleaf(tree::AbstractSimpleTree{N}, node::N) where N -> Bool

Judge whether a tree's node is a leaf (a node without children) or not.
"""
@inline isleaf(tree::AbstractSimpleTree{N}, node::N) where {N} = length(children(tree, node)) == 0

"""
    level(tree::AbstractSimpleTree{N}, node::N) where N -> Int

Get the level of tree's node.
"""
@inline function level(tree::AbstractSimpleTree{N}, node::N) where N
    result = 1
    while node != root(tree)
        result += 1
        node = parent(tree, node)
    end
    return result
end

"""
    ancestor(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N -> N

Get the ancestor of a tree's node of the n-th generation.
"""
@inline function ancestor(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N
    @assert generation >= 0 "ancestor error: generation($generation) must be non-negative."
    result = node
    for i = 1:generation
        result = parent(tree, result)
    end
    result
end

"""
    descendants(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N -> Vector{N}

Get the descendants of a tree's node of the nth generation.
"""
function descendants(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N
    @assert generation >= 0 "descendants error: generation($generation) must be non-negative."
    result = N[node]
    for i = 1:generation
        result = vcat((children(tree, node) for node in result)...)
    end
    result
end

"""
    siblings(tree::AbstractSimpleTree{N}, node::N) where N -> Vector{N}

Get the siblings (other nodes sharing the same parent) of a tree's node.
"""
@inline siblings(tree::AbstractSimpleTree{N}, node::N) where N = filter(!=(node), children(tree, parent(tree, node)))

"""
    leaves(tree::AbstractSimpleTree) -> Vector{keytype(tree)}

Get a tree's leaves.
"""
@inline leaves(tree::AbstractSimpleTree) = keytype(tree)[node for node in keys(tree, simpletreedepth) if isleaf(tree, node)]

"""
    push!(tree::AbstractSimpleTree{N, D}, node::N, data::D) where {N, D} -> typeof(tree)
    push!(tree::AbstractSimpleTree{N, D}, parent::Union{N, Nothing}, node::N, data::D) where {N, D} -> typeof(tree)

Push a new node to a tree. When `parent` is `nothing`, this function set the root node of an empty tree.
"""
Base.push!(tree::AbstractSimpleTree{N, D}, node::N, data::D) where {N, D} = push!(tree, nothing, node, data)
function Base.push!(tree::AbstractSimpleTree{N, D}, parent::Union{N, Nothing}, node::N, data::D) where {N, D}
    addnode!(tree, parent, node)
    tree[node] = data
    return tree
end

"""
    append!(tree::AbstractSimpleTree{N, D}, subtree::AbstractSimpleTree{N, D}) where {N, D} -> typeof(tree)
    append!(tree::AbstractSimpleTree{N, D}, node::Union{N, Nothing}, subtree::AbstractSimpleTree{N, D}) where {N, D} -> typeof(tree)

Append a subtree to a tree.
"""
Base.append!(tree::AbstractSimpleTree{N, D}, subtree::AbstractSimpleTree{N, D}) where {N, D} = append!(tree, nothing, subtree)
function Base.append!(tree::AbstractSimpleTree{N, D}, node::Union{N, Nothing}, subtree::AbstractSimpleTree{N, D}) where {N, D}
    for (key, value) in pairs(subtree, simpletreewidth)
        push!(tree, parent(subtree, key, node), key, value)
    end
    return tree
end

"""
    delete!(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Delete a node and all its descendants from a tree.
"""
function Base.delete!(tree::AbstractSimpleTree{N}, node::N) where N
    for key in collect(N, keys(tree, simpletreedepth, node))
        deletenode!(tree, key)
    end
    return tree
end

"""
    empty!(tree::AbstractSimpleTree) -> typeof(tree)

Empty a tree.
"""
@inline Base.empty!(tree::AbstractSimpleTree) = delete!(tree, root(tree))

"""
    subtree(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Get a subtree whose root is `node`.
"""
function subtree(tree::AbstractSimpleTree{N}, node::N) where N
    result = empty(tree)
    for (i, (key, value)) in enumerate(pairs(tree, simpletreedepth, node))
        push!(result, (i == 1) ? nothing : parent(tree, key), key, value)
    end
    return result
end

"""
    move!(tree::AbstractSimpleTree{N}, node::N, parent::N) where N -> typeof(tree)

Move a subtree to a new position.
"""
function move!(tree::AbstractSimpleTree{N}, node::N, parent::N) where N
    sub = subtree(tree, node)
    delete!(tree, node)
    append!(tree, parent, sub)
    return tree
end

"""
    SimpleTree{N, D}() where {N, D}

The minimum tree structure that implements all the default tree methods.
"""
struct SimpleTree{N, D} <: AbstractSimpleTree{N, D}
    TREECORE::SimpleTreeCore{N, D}
end
SimpleTree{N, D}() where {N, D} = SimpleTree(SimpleTreeCore{N, D}())

end #module
