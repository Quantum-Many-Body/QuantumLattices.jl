```@meta
CurrentModule = QuantumLattices
```

```@setup algorithm
using QuantumLattices
using StaticArrays: SVector, @SMatrix
ENV["GKSwstype"] = "100"  # headless output device for the GR backend of Plots
using Plots
include(joinpath(@__DIR__, "ToyTBA.jl"))
```

# [7. The Algorithm Interface](@id TutorialAlgorithmInterface)

[Chapter 6](@ref TutorialLatticeModel) ended with a [`LatticeModel`](@ref), which is a complete **description** of a quantum lattice system. Turning that description into numbers is the job of **algorithm packages** (e.g., [TightBindingApproximation](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl), [ExactDiagonalization](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl), [QuantumClusterTheories](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl), [DynamicalCorrelators](https://github.com/ZongYongyue/DynamicalCorrelators.jl), [SpinWaveTheory](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl), [MeanFieldTheory](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl), [RandomPhaseApproximation](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl), etc.), which are developed separately from [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl). The algorithm interface described in this chapter is what connects the two, and it sets itself two goals: a **uniform interface** to all the supported algorithms, so that a model is defined once and can then be handed to any of them, and **automatic project management**, so that the bookkeeping of a calculation is carried by the interface rather than left to the user.

## 7.1 Overview

### Running a minimal algorithm package

The user side of the interface is small enough to be shown in full, and it will be exercised live throughout this chapter. The example is graphene, computed with `ToyTBA`, a minimal but complete algorithm package whose whole source fits in a single file, `ToyTBA.jl`, included in this documentation and dissected piece by piece in the rest of this chapter. The real packages of Section 7.9, such as [TightBindingApproximation.jl](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl), are used in exactly the same way from the user's side.

First, build the model in the standard way of [Chapter 6](@ref TutorialLatticeModel), wrap it in the frontend, and compute the band structure along a high-symmetry path:

```@example algorithm
unitcell = Lattice([0.0, 0.0], [0.0, √3/3]; vectors=[[1.0, 0.0], [0.5, √3/2]])
hilbert = Hilbert(Fock{:f}(1, 1), length(unitcell))
model = LatticeModel(unitcell, hilbert, Hopping(:t₁, -1.0, 1))    # OperatorGenerator

graphene = Algorithm(:Graphene, ToyTBA(model))
path = ReciprocalPath(unitcell, hexagon"Γ-K-M-Γ", length=100)
energybands = graphene(:EB, EigenSystem(path))
plot(path, energybands.data.values)
```

The same computation works just as well when the Hamiltonian is supplied analytically, as a [`Formula`](@ref):

```@example algorithm
function A₀(t₁, k=SVector(0.0, 0.0); kwargs...)
    α = k[1]/2 + √3*k[2]/2
    β = -k[1]/2 + √3*k[2]/2
    h = t₁*(1 + exp(1im*α) + exp(1im*β))
    return @SMatrix [0 conj(h); h 0]
end
formula = LatticeModel(A₀, (t₁=-1.0,))                           # Formula

analytic = Algorithm(:GrapheneAnalytic, ToyTBA(formula))
energybands′ = analytic(:EB′, EigenSystem(path))
plot!(path, energybands′.data.values)
```

```@example algorithm
energybands.data.values ≈ energybands′.data.values
```

The two computations differ only in how the model is represented: an [`OperatorGenerator`](@ref) in the first, a [`Formula`](@ref) in the second. Everything else, the action `EigenSystem`, the assignments named `:EB` and `:EB′`, and the energy bands that come back, is identical. It is the frontend `ToyTBA` that absorbs the difference between the representations, so that a single algorithm works with any model.

### The five roles

Why are new structures needed at all, beyond the three representations of a [`LatticeModel`](@ref)? Consider what has to happen when a model is handed to an algorithm and a result comes back:

- The model must be recast into the form the algorithm needs, whatever representation it happens to use. Tight-binding needs a momentum-space matrix, exact diagonalization needs a many-body basis, and neither a bare [`OperatorGenerator`](@ref) nor a bare [`Formula`](@ref) knows about either. This recasting, uniform across representations, is the [`Frontend`](@ref).
- What to compute from the prepared form must be specified: an [`Action`](@ref).
- The result must have a defined shape, so that it can be stored, plotted and compared: a [`Data`](@ref).
- The method itself, with its own name and parameters, must be an object that can be updated and recorded: an [`Algorithm`](@ref).
- Each concrete task, with its own parameter values, its dependencies on other tasks, and a slot for its result, must be tracked individually: an [`Assignment`](@ref).

The roles are the same whether the algorithm diagonalizes a small Hamiltonian exactly, computes a band structure, or performs a tensor-network sweep:

* [`Frontend`](@ref): a model prepared in the form a particular algorithm needs.
* [`Action`](@ref): what to compute.
* [`Data`](@ref): the result of the computation.
* [`Algorithm`](@ref): a frontend together with the algorithm-level parameters.
* [`Assignment`](@ref): a named task that pairs an action with parameters, dependencies, and eventually its data.

A computation then proceeds as follows:

```
model
  │  Frontend(model)
  ▼
Frontend ──▶ Algorithm(name, frontend, parameters, map)
                 │
                 │  alg(name, action, [parameters,] dependencies...; options...)
                 │  ├── creates an Assignment
                 │  │     (dependencies: previously created assignments)
                 │  └── runs run!(alg, assignment), unless delay = true
                 ▼
             Assignment, now holding its Data
```

Note that [`run!`](@ref) is not a step the user performs. It is the hook an algorithm package implements, and it is invoked from inside the algorithm call shown above. A call made with `delay=true` is the exception: it returns the [`Assignment`](@ref) without computing anything, and the computation is then triggered by calling the algorithm with that assignment.

### Automatic project management

A calculation is rarely a single run: parameters are scanned, results are compared across runs, and the same model is revisited many times, often long after it was first defined. The interface therefore carries the bookkeeping along with it, rather than leaving it to the user. Every result is recorded together with the parameters that produced it, so that nothing is silently overwritten and nothing is read back wrongly. Results can be cached and reused, so that nothing is recomputed unnecessarily. And the tasks of a calculation may depend on one another, so that when the parameters of a task are changed, its next call recomputes it together with exactly those dependencies that are affected, and no others.

The two goals are precisely the two purposes that grouped the generic interface of [Section 6.1](@ref TutorialLatticeModel): the first is served by the Hamiltonian-access, parameter-management and value-type groups, and the second by the identity and persistence group, whose detailed introduction was deferred to this chapter and is given in Section 7.7.

This chapter addresses two readers at once: the user who wants to hand a model to an existing algorithm and run it, and the developer who wants to write such an algorithm.

Throughout the rest of this chapter we dissect `ToyTBA`, modelled on the real [TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) package but stripped to its essentials. It computes the band structure and the density of states of a free-fermion lattice model. By the end of the chapter every part of the protocol will have been exercised on it.

One last remark before we begin. [`Frontend`](@ref), [`Algorithm`](@ref) and [`Assignment`](@ref) are all subtypes of [`LatticeModel`](@ref). They are not further representations of a quantum lattice system in the sense of [Chapter 6](@ref TutorialLatticeModel), since none of them describes a Hamiltonian, but the common interface established there, [`str`](@ref), [`Parameters`](@ref) and [`update!`](@ref), applies to them unchanged, and the identity and file machinery built on it (Section 7.7) comes along for free. This is also what makes the second goal attainable, since the bookkeeping and the persistence are inherited by an algorithm rather than implemented by it.

## 7.2 Frontend: One Algorithm, Any Model

Let us begin with the piece that gives the algorithm its name. `ToyTBA` computes the band structure of a free-fermion lattice model, and the only thing it needs from the model is a momentum-space matrix. It should be able to get that matrix **however the model happens to be expressed**: from a [`Formula`](@ref), which returns the matrix directly, or from an [`OperatorGenerator`](@ref), whose Hamiltonian must first be turned into a matrix. (The third representation of [Chapter 6](@ref TutorialLatticeModel), [`StaticGenerator`](@ref), can be handled much like an [`OperatorGenerator`](@ref) once a [`Table`](@ref) is supplied for it; it is omitted here only for brevity.)

Neither of those two representations can do the job by itself. A [`Formula`](@ref) knows nothing about the single-particle basis, and an [`OperatorGenerator`](@ref) knows nothing about momentum space. What is missing is precisely the knowledge that belongs to the *algorithm* rather than to the *model*: which basis to use, and how a real-space operator becomes a matrix element at a given ``\mathbf{k}``. The [`Frontend`](@ref) is where that knowledge lives. It is the uniform encapsulation of one algorithm's requirements, and because it is defined for **any** [`LatticeModel`](@ref), an algorithm written against it works with every representation of a model alike.

The frontend of `ToyTBA`, quoted from `ToyTBA.jl`, is:

```julia
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
```

Note that both constructors take a [`LatticeModel`](@ref) as their only argument: the frontend is a layer above the model, whatever its representation. (Real packages often also provide constructors from the ingredients of a model, such as `TBA(lattice, hilbert, terms)`; see Section 7.9 for why.) The constructor for an [`OperatorGenerator`](@ref) reaches the Hilbert space through the `hilbert` field of the model, which stores the expanded operators, the bonds, the Hilbert space and the terms. Below the constructors, [`Base.show`](@ref) gives the frontend its label, while [`Parameters`](@ref) and [`update!`](@ref) forward to the model. The last line declares which part of the frontend enters its structural fingerprint: [`contenttoconfig`](@ref) is forwarded to the model as well, so that two `ToyTBA`s built on different models never share a configuration. The fingerprint machinery itself is the subject of Section 7.7.

The preparation required depends on the representation, which is why the two constructors differ. An analytic formula needs nothing. An [`OperatorGenerator`](@ref) does need the table, because its Hamiltonian must be assembled into a matrix by placing every operator at the row and column that its generators occupy in the single-particle basis. The real TBA package performs the same construction, and here the whole pipeline is compressed into a single function, again quoted from `ToyTBA.jl`:

```julia
matrix(frontend::ToyTBA{<:Formula}, k) = frontend.model(k)

function matrix(frontend::ToyTBA{<:OperatorGenerator}, k)
    n = length(frontend.table)
    m = zeros(ComplexF64, n, n)
    for operator in expand(frontend.model)
        m[frontend.table[operator[1]'], frontend.table[operator[2]]] += operator.value*exp(1im*dot(k, icoordinate(operator)))
    end
    return m
end
```

For each operator, the row is the basis state in which a particle is created, the column the one from which it is annihilated, and the phase is picked up from the coordinate of the operator. Here [`icoordinate`](@ref) is applied to an operator rather than to a point or a bond: for a rank-2 operator built from [`CoordinatedIndex`](@ref)es, it returns the difference of the integral coordinates of its two generators, mirroring the convention for bonds (Section 2.3.1), which is exactly the displacement that enters the Bloch phase factor. Note also that [`matrix`](@ref) is the same generic function that returned the local matrix of a spin operator in [Chapter 3](@ref TutorialDOF); the frontend simply adds its own methods to it. In this toy we only ever meet rank-2 operators on a single-orbital Fock space, so no case distinction is needed.

The graphene model and its analytic counterpart were already built in Section 7.1. Its unitcell contains two points, so the resulting matrix is a genuine ``2\times2``. Expanding the model produces 6 operators, two on each of the three 1st-neighbor bonds, and each of them carries the coordinate of the unitcell it belongs to. Those coordinates are exactly what becomes the phase factor of a Bloch state, which is why the second route of `matrix` can turn the operators into matrix elements. The two routes then give the same matrix:

```@example algorithm
numerical = ToyTBA(model)
k = SVector(pi/3, pi/2)
matrix(ToyTBA(formula), k)
```

```@example algorithm
matrix(numerical, k)
```

Note that the two methods return **different types**, an `SMatrix{2, 2}` and a `Matrix{ComplexF64}`. That difference is intentional and harmless here, because nothing in the protocol below depends on the concrete type of the matrix, only on the fact that a frontend can produce one. It is a first hint of the property we will confirm in Section 7.4, namely that an algorithm written once can work with every representation of a [`LatticeModel`](@ref).

## 7.3 Algorithm and the `map`

At this point `ToyTBA` can turn a model into a matrix, but it has no identity and no parameters of its own: `Parameters(frontend)` merely forwards to the model. The [`Algorithm`](@ref) supplies both. It wraps a frontend, gives it a name, and carries the parameters of the *method* rather than of the model. The graphene algorithm was already created in Section 7.1, and [`str`](@ref), the compact tag announced in Section 6.1, gives its label:

```@example algorithm
str(graphene)
```

With neither of the last two arguments given, the algorithm inherits the parameters of the frontend, and the translation between its own parameters and those of the frontend is the identity, which is exactly right when the two sets coincide, as they do for graphene. An algorithm is a [`LatticeModel`](@ref) in its own right, so it can be updated and inspected like any other:

```@example algorithm
update!(graphene; t₁=-0.8)
Parameters(graphene), Parameters(graphene.frontend)
```

The update acts on the frontend, and through it on the very model the frontend was built from, just as [Chapter 6](@ref TutorialLatticeModel) described for a [`LatticeModel`](@ref) and its terms.

The two parameter sets do not always coincide, however, and the [Haldane model](https://doi.org/10.1103/PhysRevLett.61.2015) is the classic counterexample. Its 2nd-neighbor hopping is complex, ``t'e^{i\varphi\pi}``, with the sign of the phase alternating between the two sublattices. A user of the algorithm thinks in terms of a magnitude ``t'`` and a flux ``\varphi``; the terms of the Hamiltonian, by contrast, need the real and imaginary parts of the hopping separately, because only the imaginary one carries the direction-dependent sign. The model therefore gains two terms, one for each part:

```@example algorithm
t₁ = Hopping(:t₁, -1.0, 1)                 # a fresh term: the one in Section 7.1 was updated to -0.8 in the update! demo
t₂ = Hopping(:t₂, -0.6*cos(0.4*pi), 2)     # the real part, with t′ = -0.6 and φ = 0.4
λ₂ = Hopping(:λ₂, -0.6*sin(0.4*pi), 2;     # the imaginary part
    amplitude=bond::Bond->1im*cos(3*azimuth(rcoordinate(bond)))*(-1)^(bond[1].site%2)
)
model = OperatorGenerator(unitcell, hilbert, (t₁, t₂, λ₂))
nothing # hide
```

In the amplitude of `λ₂`, the azimuth of a 2nd-neighbor bond of the honeycomb lattice is always a multiple of 60°, so `cos(3*azimuth(...))` evaluates to ±1 and picks out the bond direction, while `(-1)^(bond[1].site%2)` flips the sign between the two sublattices. Expanding the enriched model produces 18 operators, six on the three 1st-neighbor bonds and twelve on the six 2nd-neighbor ones; the two 2nd-neighbor terms share their generators, so the expansion combines each matching pair into a single complex operator. The analytic formula likewise gains the 2nd-neighbor contributions: `ε`, which acts equally on the two sublattices, and `d`, which acts oppositely and therefore distinguishes them.

```@example algorithm
function A(t₁, t₂, λ₂, k=SVector(0.0, 0.0); kwargs...)
    α = k[1]/2 + √3*k[2]/2
    β = -k[1]/2 + √3*k[2]/2
    h = t₁*(1 + exp(1im*α) + exp(1im*β))
    ε = 2t₂*(cos(α) + cos(β) + cos(k[1]))
    d = -2λ₂*(sin(α) - sin(β) - sin(k[1]))
    return @SMatrix [ε+d conj(h); h ε-d]
end
analytic = ToyTBA(Formula(A, (t₁=-1.0, t₂=-0.6*cos(0.4*pi), λ₂=-0.6*sin(0.4*pi))))
numerical = ToyTBA(model)
nothing # hide
```

and the two routes still give the same matrix:

```@example algorithm
matrix(analytic, k) ≈ matrix(numerical, k)
```

The algorithm, finally, is defined in the parameters of the user, and the translation into the parameters of the terms is supplied as the fourth argument, the `map`:

```@example algorithm
params(parameters::Parameters) = (t₁=parameters.t, t₂=parameters.t′*cos(parameters.φ*pi), λ₂=parameters.t′*sin(parameters.φ*pi))

algorithm = Algorithm(:Haldane, numerical, (t=-1.0, t′=-0.6, φ=0.4), params)
str(algorithm)
```

Splitting a magnitude and a phase into a real and an imaginary part is precisely what the `map` exists for. Note that the constructor does not apply the `map` by itself: the frontend is assumed to be built with the corresponding values, and from then on [`update!`](@ref) keeps the two sets in sync:

```@example algorithm
update!(algorithm; φ=0.1)
Parameters(algorithm), Parameters(algorithm.frontend)
```

Lowering the flux has moved the 2nd-neighbor hopping towards its real part, which is what a change of the enclosed flux should do.

An [`Algorithm`](@ref) is a **concrete** type, unlike [`Frontend`](@ref). Every algorithm shares the same state (a directory, a name, a frontend, a set of parameters, a `map` and a timer), so the differences between algorithms can be absorbed into the type of the frontend it holds. The reasoning behind this split is the subject of Section 7.8.

Calling the algorithm is what creates an [`Assignment`](@ref), the fifth and last role. An assignment is a *task*: it fixes one action, one set of parameter values, the dependencies it needs, and the place where its result will be put. Like the algorithm, it is a [`LatticeModel`](@ref), so it can be named, stored and plotted. We will meet it properly in Section 7.5; for now it is enough to know that it is what a call to the algorithm produces.

## 7.4 Action, Data and `run!`

A frontend produces a momentum-space matrix; what remains is to say *what to compute* from it. That is an [`Action`](@ref), and for the band structure it carries nothing but the reciprocal space to be sampled, quoted from `ToyTBA.jl`:

```julia
struct EigenSystem{R<:ReciprocalSpace} <: Action
    reciprocalspace::R
end
function Base.show(io::IO, eigensystem::EigenSystem)
    reciprocalspace = eigensystem.reciprocalspace
    reciprocalspace isa BrillouinZone && return print(io, "EigenSystem(", join(periods(reciprocalspace), "×"), ")")
    return print(io, "EigenSystem")
end
```

The reciprocal space is any subtype of [`ReciprocalSpace`](@ref): a [`BrillouinZone`](@ref) grid, as used below for the density of states, or a [`ReciprocalPath`](@ref), as used for the band-structure plots of Section 7.1. The result of the action is a [`Data`](@ref), which for the band structure is the collection of eigenvalues on that space, stored as a matrix whose rows are the sampled ``\mathbf{k}``-points and whose columns are the bands:

```julia
struct EigenSystemData <: Data
    values::Matrix{Float64}
end
```

What ties an action to its data is a single method of [`run!`](@ref), dispatched on the pair of an [`Algorithm`](@ref) and an [`Assignment`](@ref):

```julia
function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:EigenSystem}; options...)
    get(options, :showinfo, false) && @info string(assignment)
    bands = Vector{Float64}[]
    for k in assignment.action.reciprocalspace
        push!(bands, eigen(Hermitian(matrix(algorithm.frontend, k))).values)
    end
    return EigenSystemData(permutedims(reduce(hcat, bands)))
end
```

The method reads what to compute from `assignment.action`; the remaining fields of an [`Assignment`](@ref), namely its parameters, its dependencies and the slot for its data, appear in Section 7.5. The `run!` above is written once, for `Algorithm{<:ToyTBA}`, and yet the frontend it acts on may be either of the two built so far: the method reaches the matrix through `matrix(algorithm.frontend, k)` alone, without ever asking which representation of the model is behind it. The framework, for its part, never asks the developer to declare the result type by hand. Given an action and a frontend, it asks the compiler what `run!` returns:

```@example algorithm
datatype(EigenSystem, typeof(analytic))
```

```@example algorithm
datatype(EigenSystem, typeof(numerical))
```

Both give `EigenSystemData`, which is the property announced in Section 7.2: one algorithm, one `run!`, one result type, two representations of the model. Two practical consequences follow. First, `run!` must be **type-stable**, otherwise [`datatype`](@ref) cannot infer `Data` and raises an error instead of guessing; moreover, since the inference happens when an assignment is *created*, the error is raised even for a `delay=true` call that computes nothing. Second, this is why [`Data`](@ref) is an abstract type with a prescribed shape rather than a free-form container: it is the contract through which the core package, which knows nothing about any particular algorithm, can still handle its results. In particular, the fields of a `Data` are interpreted in order as the columns of a data file and as the arguments of a plot, so their order matters.

## 7.5 Creating and Running Assignments

An [`Assignment`](@ref) comes into being when the algorithm is called with an action. The call records the task under a name, and it is itself a [`LatticeModel`](@ref), so it has a directory, an identity and a place to keep its data:

```@example algorithm
eigensystem = algorithm(:eigensystem, EigenSystem(BrillouinZone(unitcell, 100)); delay=true)
```

Here the [`BrillouinZone`](@ref) is constructed directly from the lattice, exactly as taught in [Chapter 2](@ref TutorialLattice). This call was made with `delay=true`, so nothing has been computed yet. Calling the algorithm again, now with the assignment, performs the computation:

```@example algorithm
algorithm(eigensystem)
size(eigensystem.data.values)
```

The same call can also be written in the reversed order, `eigensystem(algorithm)`. The two forms differ in whose parameters take effect when the two sides disagree: `alg(assign)` first synchronizes the algorithm, and with it the frontend and the model, to the parameters of the assignment, while `assign(alg)` synchronizes the assignment to those of the algorithm. We use the first form throughout.

Assignments can also depend on one another. The density of states, for instance, is a sum over the eigenvalues computed above, so it takes that assignment as a dependency. Its action, data and `run!` are quoted from `ToyTBA.jl`:

```julia
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
```

The dependency is passed as one more positional argument of the algorithm call:

```@example algorithm
dos = algorithm(:DOS, DensityOfStates(), eigensystem)
```

and the `run!` above reads it back from `assignment.dependencies`, so the name used at the call site is never needed inside the implementation. No parameters were given to this call: an assignment inherits the current parameters of the algorithm (formally, its parameters are `merge(alg.parameters, parameters)`), and per-call values only need to be stated when the task should differ from the algorithm. (A `map` can likewise be given between the parameters and the dependencies; it defaults to `identity` and translates the parameters of the assignment into updates of its action, which none of our actions need.) This call was made without `delay`, so both the density of states and, if necessary, the eigensystem it depends on have been computed. The result is

```@example algorithm
dos.data.energies[1:3]
```

```@example algorithm
dos.data.values[1:3]
```

Since no energy range was given, it was taken from the eigenvalues of the dependency. Note also the normalization convention: every eigenvalue contributes a Gaussian of unit area, so the density of states integrates, up to the tails cut off by the energy window, to the total number of sampled single-particle states, here ``2\times100^2``; divide by the number of sampled ``\mathbf{k}`` points if the per-unitcell convention is preferred. Finally, a dependency such as this one is didactic rather than typical: the density of states of the real [TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) package computes the eigenvalues internally and takes no dependency at all.

## 7.6 Options: Self-Documenting Keyword Arguments

The algorithm call accepts keyword arguments, and an algorithm package declares which ones each of its actions understands. The declarations of `ToyTBA` are quoted from `ToyTBA.jl`:

```julia
@inline options(::Type{<:Assignment{<:EigenSystem}}) = (
    showinfo = "show the information",
)
@inline options(::Type{<:Assignment{<:DensityOfStates}}) = (
    emin = "lower bound of the energy range",
    emax = "upper bound of the energy range",
    ne = "number of sample points in the energy range",
    σ = "broadening factor",
)
```

The purpose of the declaration is to make the interface self-describing, and to let the framework police it. A call such as `alg(assign; kwargs...)` is checked against the declared options before anything else happens, so a misspelled keyword is reported together with the list of what was available:

```@example algorithm
try
    algorithm(dos; emni=-5.0)
catch error
    print(first(error.msg, 55), "...")
end
```

A user who wants the same information before writing the call can simply ask for it, and the options of the dependencies are gathered along the way:

```@example algorithm
print(optionsinfo(typeof(dos)))
```

The options declared above are not decoration. A `run!` receives them as its own keyword arguments, so `showinfo` is honoured by the `run!` of Section 7.4, and the energy range is used by the `run!` of Section 7.5. Two framework rules complete the picture. First, the check can be bypassed by passing `false` as the second positional argument of a call, which is exactly what the framework does when it runs the dependencies of a task, so the keyword arguments of a call are in fact forwarded to the `run!` of every assignment in the dependency graph. That is why the `run!` of Section 7.4 ends with `options...`: a `run!` must tolerate keywords that are meant for another action. Second, [`update!`](@ref) works on an assignment as well: it replaces the parameters of the task (pushing them through its `map` into the action), and the staleness is then detected lazily, by an approximate comparison of the parameter values at the next call. Changing a parameter therefore invalidates the task, so that the next call recomputes both it and, through the algorithm, its dependency, with the requested energy range applying to the density of states alone:

```@example algorithm
update!(dos; φ=0.3)
algorithm(dos; emin=-4.0, emax=4.0)
dos.data.energies[1:3]
```

## 7.7 Identity and Persistence: Automatic Project Management

Because [`Frontend`](@ref), [`Algorithm`](@ref) and [`Assignment`](@ref) are all [`LatticeModel`](@ref)s, the identity and persistence interfaces announced in Section 6.1 apply to them unchanged. This is not a convenience but the point of the design: it is what lets results be recorded, reused and invalidated correctly no matter which algorithm produced them. The group is introduced here, at the level where its motivation is concrete, in four parts: the in-session reuse of results, the identity of a model, the on-disk records, and the export of results.

### Not recomputing what is already known

A task remembers its parameters and its result, so a call that changes nothing costs nothing:

```@example algorithm
snapshot = eigensystem.data
algorithm(eigensystem)
eigensystem.data === snapshot
```

The identical object comes back, meaning no computation took place. Now change a parameter of the *density of states* and run it again. Its own result must be recomputed, and so must the eigensystem it depends on, even though nothing was said about the eigensystem at the call site:

```@example algorithm
update!(dos; φ=0.2)
algorithm(dos)
eigensystem.data === snapshot
```

The eigensystem data is a different object now: the dependency graph was walked, and the stale result was refreshed automatically. Two qualifications complete the picture. The walk is *top-down*: a dependency is visited only when the called task itself is stale, and it then runs with the parameters of the calling algorithm, so a calculation is steered from the top-level task, and updating a dependency directly followed by a call to a downstream task recomputes nothing. And the comparison of the parameters is approximate, so a change below the tolerance does not trigger a recomputation.

### Identity: `str`, `config` and `stamp`

[`str`](@ref) returns a compact string tag of any [`LatticeModel`](@ref), an algorithm or an assignment included; it is what labels the records and the figures below. Beyond the tag, the identity of a model has two parts: its defining structure, which is independent of the parameter values, and the parameter values themselves. The two parts exist because they serve two different distinctions: the structure is what defines the physical model, and a model remains the same model when its parameters are tuned, so records from different models must be distinguishable, and so must records from the same model at different parameter values.

The structural part is extracted by [`contenttoconfig`](@ref): for an [`OperatorGenerator`](@ref), it is the expanded coupling structure of the terms on the bonds, and for our frontend it is forwarded to the model, as quoted in Section 7.2. [`config`](@ref) returns the fingerprint of this content as a SHA-512 digest, so it is blind to parameter values, and serves the first distinction:

```@example algorithm
config(algorithm)
```

[`stamp`](@ref) appends a human-readable rendering of the current parameter values, in the form `"config-name₁(value₁)name₂(value₂)…"`, and serves the second distinction:

```@example algorithm
stamp(algorithm)
```

Therefore, [`update!`](@ref) never changes the fingerprint, and each parameter set has its own stamp. An [`Assignment`](@ref) declares no structural content of its own, so its stamp encodes its parameter values alone, within the file named after the task, which is exactly the granularity that recording needs.

### Data and cache files

The identity of a model makes its records recognizable; persistence makes them survive beyond the current session. On disk, the records live in two files, one for each kind of content: the *data* file stores the model itself, or the assignment together with its result, while the *cache* file stores the cacheable content [`contenttocache`](@ref), intermediate results that can in principle be recomputed but are expensive enough to be worth keeping (see the note at the end of this section). Both files are located in the same way: the file of a model is `pathof(model, target)`, which is `joinpath(dirname(model), basename(model, target))`. By default, [`dirname`](@ref) is `"."` (a model subtype with a `dir` field uses that directory instead), and [`basename`](@ref) is the representation name plus the extension `.qld` for `target=:data` or `.qlc` for `target=:cache`:

```@example algorithm
basename(dos, :data), basename(dos, :cache)
```

Note that the file name contains no parameter part: all parameter sets of a task share one file, and even different tasks sharing the same name coexist in it, with the stamps distinguishing the records. Within such a file, the save and load functions operate on the records:

* [`qlsave`](@ref)`(target, model)` saves a record under the key `stamp(model)`: with `target=:data` it saves the model itself into the data file, and with `target=:cache` it saves the cacheable content into the cache file. Repeated saves accumulate in the file as long as the stamps differ, while saving with an unchanged stamp replaces the old record. [`qldsave`](@ref)`(model)` and [`qlcsave`](@ref)`(model)` are the shortcuts for the two targets.
* [`qlload`](@ref) is the inverse: `qlload(path)` returns all records of a file as a `Dict` keyed by the stamps, and `qlload(path, stamp)` returns the single record under the given stamp. Therefore, an old parameter set can be recovered as long as its stamp is kept.
* [`qlclean`](@ref) removes records from a file: those not touched for longer than `age` seconds, or all but the `maxcount` most recent ones. [`qldclean`](@ref) and [`qlcclean`](@ref) apply it to all `.qld` and `.qlc` files in a directory, respectively. Every save and load refreshes the timestamp of a record, so recently used records survive the cleaning.

Saving therefore never overwrites an earlier result obtained with different parameters, and loading can never pick up the wrong one:

```@example algorithm
qldsave(dos)
oldstamp = stamp(dos)
update!(dos; φ=0.1)
algorithm(dos)
qldsave(dos)
records = qlload(pathof(dos, :data))
length(records)
```

Two records now sit in the same file, and each is retrieved by its own stamp:

```@example algorithm
qlload(pathof(dos, :data), stamp(dos)) == dos
```

```@example algorithm
Parameters(qlload(pathof(dos, :data), oldstamp))
```

The other record still carries the parameters it was computed with, untouched by the one saved afterwards. Records accumulate over time, so there is a cleanup interface as well:

```@example algorithm
qlclean(pathof(dos, :data); maxcount=1)
```

### Exporting and visualizing

A result can also be written out as plain text, in which case the fields of its `Data` are used as the columns:

```@example algorithm
dlmsave(dos)
isfile(joinpath(dirname(dos), string(str(dos), ".dlm")))
```

or visualized directly. The package ships recipes for both [Plots](https://github.com/JuliaPlots/Plots.jl) and [Makie](https://github.com/MakieOrg/Makie.jl), so with either backend loaded the same call works; and since the recipe is generic over [`Data`](@ref), it also takes care of labelling the figure with the [`str`](@ref) of the assignment:

```@example algorithm
plot(dos)
```

!!! note
    A frontend can additionally declare what should be cached through [`contenttocache`](@ref), and [`qlcsave`](@ref) then stores that content under the stamp of the model. This is how an expensive preparation, such as the expansion of a large Hamiltonian, is computed once and reused across sessions. Nothing in the toy frontend is expensive enough to be worth caching, so its default empty cache is kept.

## 7.8 Why These Five Roles

Now that all five roles have been written out and exercised, we can look back and see why the protocol has the shape it does.

The starting point is a limitation of the language. In Julia only abstract types can be supertypes, and abstract types have no fields. So the state shared by every algorithm, or by every task, cannot be hoisted into a common base class; it has to live in **concrete** types. The framework therefore splits the two concerns: abstract types carry the *behaviour contract*, since methods defined on them are inherited by every subtype, while concrete types carry the *state*. This is why [`Frontend`](@ref) and [`Action`](@ref) are abstract while [`Algorithm`](@ref) and [`Assignment`](@ref) are concrete, and it is why the latter two are parameterised by the former:

| Role | Kind | State | Why |
|------|------|-------|-----|
| [`Frontend`](@ref) | abstract | not shared | every package's frontend holds something different |
| [`Action`](@ref) | abstract | not shared | every action specifies something different |
| [`Data`](@ref) | abstract | not shared | every result has a different shape |
| [`Algorithm`](@ref) | concrete | shared | every algorithm has the same directory, name, frontend, parameters, `map` and timer |
| [`Assignment`](@ref) | concrete | shared | every task has the same directory, name, action, parameters, `map`, dependencies and data |

The two concrete types divide the labour between "how to compute" and "what is being computed". An [`Algorithm`](@ref) fixes the method and the parameter values; an [`Assignment`](@ref) fixes a single task, and earns its place by three things an algorithm alone could not provide. It is **named and located**, so its result can be saved, plotted and cited. It keeps its **data in an uninitialised slot**, so that a task can be described without being run, and re-run only when its parameters have changed. And it holds its **dependencies**, which turn a collection of tasks into a computation graph that the framework walks in the right order.

That leaves the question of why [`Data`](@ref) is designed the way it is, and the answer is subtler than it first appears.

The decisive reason is that an [`Assignment`](@ref) is constructed **before** anything is computed, and may never be computed at all if `delay=true`. Its result slot must nevertheless have a concrete type, and Julia allows a field to be left undefined only if its type is known at construction time. The type of the result therefore cannot be decided after the fact; it has to be derived statically from the pair of the frontend and the action, which is exactly what [`datatype`](@ref) does by asking the compiler for the return type of [`run!`](@ref). This is why `Data` must be an abstract type with a well-defined interface, why the type of the result is a type parameter of the [`Assignment`](@ref) rather than a property of the action, and why `run!` is required to be type-stable.

The second reason is genericity of the operations that come after the computation. The `Tuple` of a `Data` returns its fields in order, and both the plain-text output and the plotting recipes are defined in terms of that tuple. The core package can therefore save and draw the results of an algorithm it has never heard of, in the same way that it can store and stamp a model it has never heard of.

The third reason is that a result is not a description. Unlike the other three types, a `Data` has no parameters, no configuration, and no identity, so there is nothing to inherit from [`LatticeModel`](@ref) and no reason to pretend otherwise. It keeps only what is needed: equality for comparison, and the tuple view for output.

## 7.9 The Ecosystem

Real algorithms are developed in separate packages, all of them built on the protocol of this chapter. The following table lists the main ones; the same packages are listed on the home page of these docs, which always reflects the current state of the ecosystem.

| Package | Method | Applicable systems |
|---------|--------|--------------------|
| [TightBindingApproximation.jl](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) | Tight-binding approximation | Fermions, bosons, phonons |
| [ExactDiagonalization.jl](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl) | Exact diagonalization | Fermions, hard-core bosons, spins |
| [QuantumClusterTheories.jl](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl) | Cluster perturbation theory, variational cluster approach | Fermions, spins |
| [DynamicalCorrelators.jl](https://github.com/ZongYongyue/DynamicalCorrelators.jl) | Density matrix renormalization group | Fermions, hard-core bosons, spins |
| [SpinWaveTheory.jl](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl) | Linear spin wave theory | Ordered spin systems |
| [MeanFieldTheory.jl](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl) | Self-consistent mean field theory | Fermions |
| [RandomPhaseApproximation.jl](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl) | Random phase approximation | Fermions |

Their frontends are richer than the one built here, and some of them deal with matrix sizes that are not known until the model is expanded, but the protocol is the same: prepare the model, say what to compute, declare the options, implement `run!`, and let the framework take care of naming, recording and reusing the results. Two pieces of the real packages are worth knowing in advance. Their frontends are usually constructed directly from the ingredients of a model, e.g., `TBA(lattice, hilbert, terms)`, which is more than a convenience: the ingredient constructor is where algorithm-specific preparation options live, such as the boundary twists of a tight-binding calculation, and such options do not belong to the model itself. `ToyTBA` omits this constructor to keep its interface uniform. And the protocol is not only uniform but *inherited*: [SpinWaveTheory.jl](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl) defines its frontend as a subtype of TBA's and reuses its actions wholesale.

## Summary

The algorithm interface connects a [`LatticeModel`](@ref) to a numerical method:

- A [`Frontend`](@ref) encapsulates what one algorithm needs from a model, and because it is written against [`LatticeModel`](@ref), a single algorithm works with every representation of a model.
- An [`Action`](@ref) and a [`Data`](@ref) say what to compute and what comes out, and the two are wired together automatically, [`datatype`](@ref) inferring the latter from the return type of `run!`.
- An [`Algorithm`](@ref) holds the method and its parameters, with a `map` translating them into the parameters of the frontend when the two sets do not coincide.
- Calling the algorithm creates an [`Assignment`](@ref), which names the task, remembers its dependencies, and holds its result; with `delay=true` the call returns the assignment without computing anything.
- [`options`](@ref) and [`optionsinfo`](@ref) make the keyword arguments of a call self-documenting and checkable.
- Because all three are [`LatticeModel`](@ref)s, [`config`](@ref) and [`stamp`](@ref) give every result a stable key, so saving never overwrites a result obtained with different parameter values and loading by stamp never confuses them, while unchanged results are never recomputed and stale dependencies are refreshed automatically.

The complete source of the toy package dissected in this chapter is the file `ToyTBA.jl` included in this documentation, and the docstrings of the types and functions introduced in this chapter are collected in the [manual](@ref "Frameworks").
