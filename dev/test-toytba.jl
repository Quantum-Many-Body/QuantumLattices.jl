# Standalone test for docs/src/tutorials/ToyTBA.jl
# Run: julia --project=docs dev/test-toytba.jl

include(joinpath(@__DIR__, "..", "docs", "src", "tutorials", "ToyTBA.jl"))
using QuantumLattices
using QuantumLattices: datatype, dependencytypes, optionsinfo
using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: eigen, Hermitian

# --- 1. graphene model, OperatorGenerator route (ch7 7.2) ---
lattice = Lattice((0.0, 0.0), (0.0, sqrt(3)/3); vectors=[[1.0, 0.0], [0.5, sqrt(3)/2]], name=:Honeycomb)
hilbert = Hilbert(Fock{:f}(1, 1), length(lattice))
t₁ = Hopping(:t₁, -1.0, 1)
model = OperatorGenerator(lattice, hilbert, t₁)
numerical = ToyTBA(model)
@assert length(expand(model)) == 6

# --- 2. analytic route (Formula) ---
function A₀(t₁, k=SVector(0.0, 0.0); kwargs...)
    α = k[1]/2 + √3*k[2]/2
    β = -k[1]/2 + √3*k[2]/2
    h = t₁*(1 + exp(1im*α) + exp(1im*β))
    return @SMatrix [0 conj(h); h 0]
end
analytic = ToyTBA(Formula(A₀, (t₁=-1.0,)))

k = SVector(pi/3, pi/2)
ma = matrix(analytic, k)
mn = matrix(numerical, k)
@assert ma ≈ mn

# --- 3. contract: valtype / dependencytypes / @delegate ---
@assert valtype(analytic) == valtype(typeof(analytic)) == SMatrix{2, 2, ComplexF64, 4}
@assert valtype(numerical) == valtype(typeof(numerical)) == Matrix{ComplexF64}
@assert dependencytypes(DensityOfStates) == (EigenSystem,)
@assert dependencytypes(EigenSystem) |> isnothing

# --- 4. Algorithm, update!, Parameters ---
algorithm = Algorithm(:Graphene, numerical)
@assert Parameters(algorithm) == (t₁=-1.0,)
update!(algorithm; t₁=-0.8)
@assert Parameters(algorithm.frontend) == (t₁=-0.8,)
@assert matrix(algorithm, k) ≈ matrix(numerical, k)   # delegated query: matrix(alg, k) ≡ matrix(alg.frontend, k)

# --- 5. datatype inference ---
@assert datatype(typeof(EigenSystem(BrillouinZone(lattice, 2))), typeof(analytic)) <: EigenSystemData
@assert datatype(typeof(EigenSystem(BrillouinZone(lattice, 2))), typeof(numerical)) <: EigenSystemData

# --- 6. assignments: delay, run, reversed call ---
eigensystem = algorithm(:eigensystem, EigenSystem(BrillouinZone(lattice, 20)); delay=true)
@assert !isdefined(eigensystem, :data)
algorithm(eigensystem)
@assert size(eigensystem.data.values) == (400, 2)

# --- 7. DOS with dependency ---
dos = algorithm(:DOS, DensityOfStates(), eigensystem)
@assert length(dos.data.energies) == 101
@assert size(dos.data.values) == (101, 1)

# --- 8. options machinery ---
@assert optionsinfo(typeof(dos)) isa String
let caught = false
    try
        algorithm(dos; emni=-5.0)
    catch
        caught = true
    end
    @assert caught
end

# --- 9. update! staleness + recompute; definitional knobs live in the task ---
snapshot = eigensystem.data
algorithm(eigensystem)
@assert eigensystem.data === snapshot
update!(dos; t₁=-0.9)
algorithm(dos)
dos₂ = algorithm(:DOS₂, DensityOfStates(emin=-4.0, emax=4.0), eigensystem)
@assert dos₂.data.energies[1] ≈ -4.0
@assert length(dos₂.data.energies) == 101

# --- 10. save/load/clean (in temp dir) ---
mktempdir() do dir
    cd(dir) do
        alg = Algorithm(:Graphene, ToyTBA(LatticeModel(lattice, hilbert, t₁)); dir=dir)
        es = alg(:eigensystem, EigenSystem(BrillouinZone(lattice, 10)))
        d = alg(:DOS, DensityOfStates(), es)
        qldsave(d)
        oldstamp = stamp(d)
        oldparams = Parameters(d)
        update!(d; t₁=-0.85)
        alg(d)
        qldsave(d)
        path = pathof(d, :data)
        @assert length(qlload(path)) == 2
        @assert isequal(qlload(path, stamp(d)), d)   # isequal, not ==: the unset energy window of the task is NaN
        @assert Parameters(qlload(path, oldstamp)) == oldparams
        qlclean(path; maxcount=1)
        @assert length(qlload(path)) == 1
    end
end

println("ALL TESTS PASSED")
