#=
    ToyTBA

A minimal but complete algorithm package built on the algorithm interface of QuantumLattices.
It computes the band structure (EigenSystem) and the density of states (DensityOfStates)
of a free-fermion lattice model.

This file is included by the documentation of QuantumLattices and dissected piece by piece
in Chapter 8 of the tutorial. In a real package, the definitions below would be wrapped
in a module; here they are meant to be `include`d directly.
=#

using QuantumLattices
import QuantumLattices: Parameters, contenttoconfig, dependencytypes, matrix, options, run!, update!
using LinearAlgebra: Hermitian, dot, eigen

# The frontend role: a model prepared in the form this algorithm needs

struct ToyTBA{M<:LatticeModel, T<:Union{Table, Nothing}} <: LatticeModel
    hamiltonian::M
    table::T
end
ToyTBA(hamiltonian::Formula) = ToyTBA(hamiltonian, nothing)
ToyTBA(hamiltonian::OperatorGenerator) = ToyTBA(hamiltonian, Table(hamiltonian.hilbert, OperatorIndexToTuple(:site, :orbital, :spin)))

@inline Base.valtype(::Type{<:ToyTBA{<:Formula{V}}}) where V = V
@inline Base.valtype(::Type{<:ToyTBA{<:OperatorGenerator}}) = Matrix{ComplexF64}
@inline Base.show(io::IO, ::ToyTBA) = print(io, "ToyTBA")
@inline Parameters(toytba::ToyTBA) = Parameters(toytba.hamiltonian)
@inline update!(toytba::ToyTBA; parameters...) = (update!(toytba.hamiltonian; parameters...); toytba)
@inline contenttoconfig(toytba::ToyTBA) = contenttoconfig(toytba.hamiltonian)

# The momentum-space matrix, one method per representation of the model

@delegate matrix(toytba::ToyTBA{<:Formula}, k) = toytba.hamiltonian(k)

@delegate function matrix(toytba::ToyTBA{<:OperatorGenerator}, k)
    n = length(toytba.table)
    m = zeros(ComplexF64, n, n)
    for operator in expand(toytba.hamiltonian)
        m[toytba.table[operator[1]'], toytba.table[operator[2]]] += operator.value*exp(1im*dot(k, icoordinate(operator)))
    end
    return m
end

# Task and Data: the band structure

struct EigenSystem{R<:ReciprocalSpace}
    reciprocalspace::R
end
function Base.show(io::IO, eigensystem::EigenSystem)
    reciprocalspace = eigensystem.reciprocalspace
    reciprocalspace isa BrillouinZone && return print(io, "EigenSystem(", join(periods(reciprocalspace), "×"), ")")
    return print(io, "EigenSystem")
end

struct EigenSystemData{R<:ReciprocalSpace} <: Data
    reciprocalspace::R
    values::Matrix{Float64}
end

function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:EigenSystem}; options...)
    get(options, :showinfo, false) && @info string(assignment)
    bands = Vector{Float64}[]
    for k in assignment.task.reciprocalspace
        push!(bands, eigen(Hermitian(matrix(algorithm, k))).values)
    end
    return EigenSystemData(assignment.task.reciprocalspace, permutedims(reduce(hcat, bands)))
end

@inline options(::Type{<:Assignment{<:EigenSystem}}) = (
    showinfo = "show the information",
)

# Task and Data: the density of states

Base.@kwdef struct DensityOfStates
    emin::Float64 = NaN
    emax::Float64 = NaN
    ne::Int = 101
    σ::Float64 = 0.1
end
@inline dependencytypes(::Type{<:DensityOfStates}) = (EigenSystem,)

struct DensityOfStatesData <: Data
    energies::Vector{Float64}
    values::Matrix{Float64}
end

function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:DensityOfStates}; options...)
    eigensystem = first(assignment.dependencies)
    (; emin, emax, ne, σ) = assignment.task
    isnan(emin) && (emin = minimum(eigensystem.data.values))
    isnan(emax) && (emax = maximum(eigensystem.data.values))
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
