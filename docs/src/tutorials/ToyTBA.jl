#=
    ToyTBA

A minimal but complete algorithm package built on the algorithm interface of QuantumLattices.
It computes the band structure (EigenSystem) and the density of states (DensityOfStates)
of a free-fermion lattice model.

This file is included by the documentation of QuantumLattices and dissected piece by piece
in Chapter 7 of the tutorial. In a real package, the definitions below would be wrapped
in a module; here they are meant to be `include`d directly.
=#

using QuantumLattices
import QuantumLattices: Parameters, contenttoconfig, matrix, options, run!, update!
using LinearAlgebra: Hermitian, dot, eigen

# Frontend: a model prepared in the form this algorithm needs

struct ToyTBA{M<:LatticeModel, T<:Union{Table, Nothing}} <: Frontend
    model::M
    table::T
end
ToyTBA(model::Formula) = ToyTBA(model, nothing)
ToyTBA(model::OperatorGenerator) = ToyTBA(model, Table(model.hilbert, OperatorIndexToTuple(:site, :orbital, :spin)))

@inline Base.show(io::IO, ::ToyTBA) = print(io, "ToyTBA")
@inline Parameters(frontend::ToyTBA) = Parameters(frontend.model)
@inline update!(frontend::ToyTBA; parameters...) = (update!(frontend.model; parameters...); frontend)
@inline contenttoconfig(frontend::ToyTBA) = contenttoconfig(frontend.model)

# The momentum-space matrix, one method per representation of the model

matrix(frontend::ToyTBA{<:Formula}, k) = frontend.model(k)

function matrix(frontend::ToyTBA{<:OperatorGenerator}, k)
    n = length(frontend.table)
    m = zeros(ComplexF64, n, n)
    for operator in expand(frontend.model)
        m[frontend.table[operator[1]'], frontend.table[operator[2]]] += operator.value*exp(1im*dot(k, icoordinate(operator)))
    end
    return m
end

# Action and Data: the band structure

struct EigenSystem{R<:ReciprocalSpace} <: Action
    reciprocalspace::R
end
function Base.show(io::IO, eigensystem::EigenSystem)
    reciprocalspace = eigensystem.reciprocalspace
    reciprocalspace isa BrillouinZone && return print(io, "EigenSystem(", join(periods(reciprocalspace), "×"), ")")
    return print(io, "EigenSystem")
end

struct EigenSystemData <: Data
    values::Matrix{Float64}
end

function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:EigenSystem}; options...)
    get(options, :showinfo, false) && @info string(assignment)
    bands = Vector{Float64}[]
    for k in assignment.action.reciprocalspace
        push!(bands, eigen(Hermitian(matrix(algorithm.frontend, k))).values)
    end
    return EigenSystemData(permutedims(reduce(hcat, bands)))
end

@inline options(::Type{<:Assignment{<:EigenSystem}}) = (
    showinfo = "show the information",
)

# Action and Data: the density of states

struct DensityOfStates <: Action end

struct DensityOfStatesData <: Data
    energies::Vector{Float64}
    values::Matrix{Float64}
end

function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:DensityOfStates}; emin=nothing, emax=nothing, ne::Int=101, σ=0.1)
    @assert isa(assignment.dependencies, Tuple{Assignment{<:EigenSystem}}) "run! error: wrong dependencies."
    eigensystem = first(assignment.dependencies)
    isnothing(emin) && (emin = minimum(eigensystem.data.values))
    isnothing(emax) && (emax = maximum(eigensystem.data.values))
    data = DensityOfStatesData(collect(range(emin, emax, ne)), zeros(ne, 1))
    for (i, ω) in enumerate(data.energies)
        data.values[i] = 0.0
        for energy in eigensystem.data.values
            data.values[i] += exp(-(ω-energy)^2/2σ^2)
        end
        data.values[i] /= √(2pi)*σ
    end
    return data
end

@inline options(::Type{<:Assignment{<:DensityOfStates}}) = (
    emin = "lower bound of the energy range",
    emax = "upper bound of the energy range",
    ne = "number of sample points in the energy range",
    σ = "broadening factor",
)
