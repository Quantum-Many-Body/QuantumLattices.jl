module Combinatorics

export AbstractCombinatorics
export Combinations, DulCombinations, Permutations, DulPermutations

"""
    AbstractCombinatorics{M, C}

Abstract combinatorial algorithms.
"""
abstract type AbstractCombinatorics{M, C} end
@inline Base.eltype(::Type{<:AbstractCombinatorics{M, C}}) where {M, C} = NTuple{M, eltype(C)}

"""
    Combinations{M}(contents::C) where {M, C}

Combinations of M elements from contents. Duplicates are not allowed.
"""
struct Combinations{M, C} <: AbstractCombinatorics{M, C}
    contents::C
    N::Int
    Combinations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(c::Combinations{M}) where M = binomial(c.N, M)
Base.iterate(c::Combinations{M}) where M = (M > c.N) ? nothing : (M == 0) ? ((), [c.N+2]) : (ntuple(i->c.contents[i], Val(M)), nextmstate!(collect(1:M), c.N, M))
Base.iterate(c::Combinations{M}, state) where M = (state[1] > c.N-M+1) ? nothing : (ntuple(i->c.contents[state[i]], Val(M)), nextmstate!(state, c.N, M))
function nextmstate!(state::Vector{Int}, N::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        (state[i] > N-(M-i)) && continue
        for j = i+1:M
            state[j] = state[j-1] + 1
        end
        break
    end
    state
end

"""
    DulCombinations{M}(contents::C) where {M, C}

Combinations of M elements from contents. Duplicates are allowed.
"""
struct DulCombinations{M, C} <: AbstractCombinatorics{M, C}
    contents::C
    N::Int
    DulCombinations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(c::DulCombinations{M}) where M = binomial(c.N+M-1, c.N-1)
function Base.iterate(c::DulCombinations{M}) where M
    (M == 0) ? ((), [c.N+1]) : (ntuple(i->c.contents[1], Val(M)), nextdmstate!(collect(ntuple(i->1, M)), c.N, M))
end
Base.iterate(c::DulCombinations{M}, state) where M = (state[1] > c.N) ? nothing : (ntuple(i->c.contents[state[i]], Val(M)), nextdmstate!(state, c.N, M))
function nextdmstate!(state::Vector{Int}, N::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        (state[i] > N) && continue
        for j = i+1:M
            state[j] = state[j-1]
        end
        break
    end
    state
end

"""
    Permutations{M}(contents::C) where {M, C}

Permutations of M elements from contents. Duplicates are not allowed.
"""
struct Permutations{M, C} <: AbstractCombinatorics{M, C}
    contents::C
    N::Int
    Permutations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(p::Permutations{M}) where M = (0 <= M <= p.N) ? prod((p.N-M+1):p.N) : 0
function Base.iterate(p::Permutations{M}) where M
    ((p.N == 0) && (M > 0) || (0 < p.N < M)) ? nothing : (state = collect(1:p.N); (ntuple(i->p.contents[state[i]], Val(M)), nextpstate!(state, p.N, M)))
end
function Base.iterate(p::Permutations{M}, state) where M
    ((p.N == 0) && (M > 0) || (0 < p.N < max(state[1], M))) ? nothing : (ntuple(i->p.contents[state[i]], Val(M)), nextpstate!(state, p.N, M))
end
function nextpstate!(state::Vector{Int}, N::Int, M::Int)
    (M <= 0) && return [N+1]
    (M < N) && (j = M + 1; while (j <= N) && (state[M] >= state[j]); j += 1; end)
    if (M < N) && (j <= N)
        state[M], state[j] = state[j], state[M]
    else
        (M < N) && reverse!(state, M+1)
        i = M - 1
        while (i >= 1) && (state[i] >= state[i+1]); i -= 1; end
        if i > 0
            j = N
            while (j > i) && (state[i] >= state[j]); j -= 1; end
            state[i], state[j] = state[j], state[i]
            reverse!(state, i+1)
        else
            state[1] = N + 1
        end
    end
    return state
end

"""
    DulPermutations{M}(contents::C) where {M, C}

Permutations of M elements from contents. Duplicates are allowed.
"""
struct DulPermutations{M, C} <: AbstractCombinatorics{M, C}
    contents::C
    N::Int
    DulPermutations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(p::DulPermutations{M}) where M = p.N^M
function Base.iterate(p::DulPermutations{M}) where M
    indices = CartesianIndices(ntuple(i->p.N, M|>Val))
    index = iterate(indices)
    isnothing(index) && return nothing
    return ntuple(i->p.contents[reverse(index[1].I)[i]], M|>Val), (indices, index[2])
end
function Base.iterate(p::DulPermutations{M}, state) where M
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return ntuple(i->p.contents[reverse(index[1].I)[i]], M|>Val), (state[1], index[2])
end

end # module
