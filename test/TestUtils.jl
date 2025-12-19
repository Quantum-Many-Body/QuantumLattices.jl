module TestUtils

using LinearAlgebra: dot
using Printf: @printf
using QuantumLattices.QuantumOperators
using QuantumLattices.Toolkit
using QuantumLattices.DegreesOfFreedom
using QuantumLattices.Spatials: decompose

import QuantumLattices: permute, shape
import QuantumLattices.DegreesOfFreedom: internalindextype, isdefinite, statistics
import QuantumLattices.QuantumOperators: latexname, script

export AID, DID, DFock, 𝕕

# AID
struct AID{O<:Real, S<:Real} <: OperatorIndex
    orbital::O
    nambu::S
end
@inline Base.adjoint(id::AID) = AID(id.orbital, 3-id.nambu)
@inline script(id::AID, ::Val{:orbital}; kwargs...) = id.orbital
@inline script(id::AID, ::Val{:nambu}; kwargs...) = id.nambu==2 ? "\\dagger" : ""
latexformat(AID, LaTeX{(:nambu,), (:orbital,)}('c'))
function permute(u₁::AID, u₂::AID)
    @assert u₁ ≠ u₂ "permute error: permuted operator units should not be equal to each other."
    if (u₁.nambu == 3-u₂.nambu) && (u₁.orbital == u₂.orbital)
        if u₁.nambu == 2
            return (Operator(1), Operator(1, u₂, u₁))
        else
            return (Operator(-1), Operator(1, u₂, u₁))
        end
    else
        return (Operator(1, u₂, u₁),)
    end
end

# DID
struct DID{N<:Union{Int, Symbol, Colon}} <: InternalIndex
    nambu::N
end

@inline Base.show(io::IO, ::Type{<:DID}) = @printf io "%s" "DID"
@inline Base.adjoint(sl::DID{Int}) = DID(3-sl.nambu)
@inline statistics(::Type{<:DID}) = :f
function permute(did₁::DID, did₂::DID)
    @assert did₁ ≠ did₂ "permute error: two identical fermionic indexes should vanish due to the fermionic statistics."
    return (Operator(1), Operator(-1, did₂, did₁))
end
@inline isdefinite(::Type{DID{Int}}) = true
@inline script(did::DID, ::Val{:nambu}; kwargs...) = did.nambu==Colon() ? ":" : string(did.nambu)
function Base.angle(id::CoordinatedIndex{<:Index{DID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float64}}, values::AbstractVector{Float64})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoordinate, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.internal.nambu == 1) ? phase : -phase
end
@inline internalindextype(::Type{DID}, ::Type{T}) where {T<:Union{Int, Symbol, Colon}} = DID{T}
@inline 𝕕(nambu) = DID(nambu)
@inline 𝕕(site, nambu) = Index(site, DID(nambu))
@inline 𝕕(site, nambu, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, DID(nambu)), rcoordinate, icoordinate)
@inline Base.getindex(::Type{OperatorIndex}, ::DID) = "𝕕"

@inline latexname(::Type{<:CoordinatedIndex{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}}) = Symbol("CoordinatedIndex{Index{DID, Union{Int, Ordinal, Colon}}}")
@inline latexname(::Type{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}) = Symbol("Index{DID, Union{Int, Ordinal, Colon}}")

latexformat(CoordinatedIndex{<:Index{<:DID, <:Union{Int, Ordinal, Colon}}}, LaTeX{(), (:site, :nambu)}('d'))
latexformat(Index{<:DID, <:Union{Int, Ordinal, Colon}}, LaTeX{(), (:site, :nambu)}('d'))
latexformat(DID, LaTeX{(), (:nambu,)}('d'))

# DFock
struct DFock <: SimpleInternal{DID{Int}}
    nnambu::Int
end

@inline shape(f::DFock) = (1:f.nnambu,)
@inline Base.show(io::IO, ::Type{DFock}) = @printf io "%s" "DFock"
@inline Base.convert(::Type{<:DID}, i::CartesianIndex, ::DFock) = DID(i.I...)
@inline Base.convert(::Type{<:CartesianIndex}, did::DID{Int}, ::DFock) = CartesianIndex(did.nambu)
@inline shape(::DFock, index::DID{Int}) = (index.nambu:index.nambu,)

end
