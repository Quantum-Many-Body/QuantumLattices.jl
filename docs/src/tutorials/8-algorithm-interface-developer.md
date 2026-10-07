```@meta
CurrentModule = QuantumLattices
```

```@setup algorithmdev
using QuantumLattices
using StaticArrays: SVector, @SMatrix
include(joinpath(@__DIR__, "ToyTBA.jl"))

unitcell = Lattice([0.0, 0.0], [0.0, √3/3]; vectors=[[1.0, 0.0], [0.5, √3/2]])
hilbert = Hilbert(Fock{:f}(1, 1), length(unitcell))
model = LatticeModel(unitcell, hilbert, Hopping(:t, -1.0, 1)) # OperatorGenerator

function A₀(t, k=SVector(0.0, 0.0); kwargs...)
    α = k[1]/2 + √3*k[2]/2
    β = -k[1]/2 + √3*k[2]/2
    h = t*(1 + exp(1im*α) + exp(1im*β))
    return @SMatrix [0 conj(h); h 0]
end
model′ = LatticeModel(A₀, (t=-1.0,)) # Formula

numerical = ToyTBA(model)
analytic = ToyTBA(model′)
path = ReciprocalPath(unitcell, hexagon"Γ-K-M-Γ", length=100)
k = SVector(pi/3, pi/2)
```

# [8. The Algorithm Interface: Developer Guide](@id TutorialAlgorithmInterfaceDeveloper)

## 8.1 Dissecting a Minimal Algorithm Package

[Chapter 7](@ref TutorialAlgorithmInterface) showed the algorithm interface from the user's side; this chapter shows it from the developer's side, by dissecting `ToyTBA`, the minimal but complete algorithm package that powered every example of that chapter. Modelled on the real [TightBindingApproximation](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) package, its whole source fits in a single file, `ToyTBA.jl`, included in this documentation, and the dissection below follows the order of that file. It computes the band structure and the density of states of a free-fermion lattice model, and by the end of the chapter every part of the protocol will have been exercised on it.

Writing an algorithm package means fulfilling a small contract:

* a **frontend**, a [`LatticeModel`](@ref) subtype that recasts any model into the form the algorithm needs, together with the query methods the algorithm will call, such as `matrix` for a momentum-space method;
* **tasks**, plain structs saying *what* is computed, and a [`Data`](@ref) subtype per result, wired together by [`run!`](@ref) methods dispatched on the pair of an [`Algorithm`](@ref) and an [`Assignment`](@ref);
* [`dependencytypes`](@ref) declarations for the tasks that depend on other assignments;
* [`options`](@ref) declarations for the keyword arguments each task understands.

Everything else is supplied by the framework: the [`Algorithm`](@ref) and [`Assignment`](@ref) types, the inference of the result type, the validation at registration time, the parameter management and the persistence. The running example is graphene again: the model and its analytic counterpart of [Chapter 7](@ref TutorialAlgorithmInterface) are rebuilt in the setup of this page, wrapped by the two frontends `numerical` and `analytic`.

## 8.2 The Frontend

The first piece of `ToyTBA.jl` is the one that gives the algorithm its name. `ToyTBA` computes the band structure of a free-fermion lattice model, and the only thing it needs from the model is a momentum-space matrix. It should be able to get that matrix **however the model happens to be expressed**: from a [`Formula`](@ref), which returns the matrix directly, or from an [`OperatorGenerator`](@ref), whose Hamiltonian must first be turned into a matrix. (The third representation of [Chapter 6](@ref TutorialLatticeModel), [`StaticGenerator`](@ref), can be handled much like an [`OperatorGenerator`](@ref) once a [`Table`](@ref) is supplied for it; it is omitted here only for brevity.)

Neither of those two representations can do the job by itself. A [`Formula`](@ref) knows nothing about the single-particle basis, and an [`OperatorGenerator`](@ref) knows nothing about momentum space. What is missing is precisely the knowledge that belongs to the *algorithm* rather than to the *model*: which basis to use, and how a real-space operator becomes a matrix element at a given ``\mathbf{k}``. The frontend is where that knowledge lives. It is the uniform encapsulation of one algorithm's requirements, and because it is defined for **any** [`LatticeModel`](@ref), an algorithm written against it works with every representation of a model alike.

The frontend of `ToyTBA`, quoted from `ToyTBA.jl`, is:

```julia
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
```

Note that both constructors take a [`LatticeModel`](@ref) as their only argument: the frontend is a layer above the model, whatever its representation. (Real packages often also provide constructors from the ingredients of a model, such as `TBA(lattice, hilbert, terms)`; see [Chapter 7](@ref TutorialAlgorithmInterface) for why.) The constructor for an [`OperatorGenerator`](@ref) reaches the Hilbert space through the `hilbert` field of the wrapped Hamiltonian, which stores the bonds, the Hilbert space and the terms. The two `valtype` lines below the constructors fulfil the one obligation that the [`LatticeModel`](@ref) supertype imposes on its subtypes: a type-level [`valtype`](@ref), here the type of the momentum-space matrix. For a [`Formula`](@ref) it is the formula's own return type, and for an [`OperatorGenerator`](@ref) it is a dense `Matrix{ComplexF64}`, since the Bloch phases make the matrix complex even when the operators are real. Then [`Base.show`](@ref) gives the frontend its label, while [`Parameters`](@ref) and [`update!`](@ref) forward to the Hamiltonian. The last line declares which part of the frontend enters its structural fingerprint: [`contenttoconfig`](@ref) is forwarded as well, so that two `ToyTBA`s built on different models never share a configuration. The fingerprint machinery itself is presented in [Chapter 7](@ref TutorialAlgorithmInterface).

The preparation required depends on the representation, which is why the two constructors differ. An analytic formula needs nothing. An [`OperatorGenerator`](@ref) does need the table, because its Hamiltonian must be assembled into a matrix by placing every operator at the row and column that its generators occupy in the single-particle basis. The real TBA package performs the same construction, and here the whole pipeline is compressed into a single function, again quoted from `ToyTBA.jl`:

```julia
@delegate matrix(toytba::ToyTBA{<:Formula}, k) = toytba.hamiltonian(k)

@delegate function matrix(toytba::ToyTBA{<:OperatorGenerator}, k)
    n = length(toytba.table)
    m = zeros(ComplexF64, n, n)
    for operator in expand(toytba.hamiltonian)
        m[toytba.table[operator[1]'], toytba.table[operator[2]]] += operator.value*exp(1im*dot(k, icoordinate(operator)))
    end
    return m
end
```

For each operator, the row is the basis state in which a particle is created, the column the one from which it is annihilated, and the phase is picked up from the coordinate of the operator. Here [`icoordinate`](@ref) is applied to an operator rather than to a point or a bond: for a rank-2 operator built from [`CoordinatedIndex`](@ref)es, it returns the difference of the integral coordinates of its two generators, mirroring the convention for bonds ([Section 2.3.1](@ref TutorialBond)), which is exactly the displacement that enters the Bloch phase factor. Note also that [`matrix`](@ref) is the same generic function that returned the local matrix of a spin operator in [Chapter 3](@ref TutorialDOF); the frontend simply adds its own methods to it. In this toy we only ever meet rank-2 operators on a single-orbital Fock space, so no case distinction is needed.

The [`@delegate`](@ref) macro marks these methods as part of the frontend's *query interface*: each of them is automatically mirrored on [`Algorithm`](@ref), so that `matrix(alg, k)` means `matrix(alg.frontend, k)`. An algorithm is a transparent proxy of its frontend for such queries; the only thing it ever injects is execution context, such as the timer.

The graphene unitcell contains two points, so the resulting matrix is a genuine ``2\times2``. Expanding the model produces 6 operators, two on each of the three 1st-neighbor bonds, and each of them carries the coordinate of the unitcell it belongs to. Those coordinates are exactly what becomes the phase factor of a Bloch state, which is why the second route of `matrix` can turn the operators into matrix elements. The two routes then give the same matrix:

```@example algorithmdev
matrix(analytic, k)
```

```@example algorithmdev
matrix(numerical, k)
```

Note that the two methods return **different types**, an `SMatrix{2, 2}` and a `Matrix{ComplexF64}`. That difference is intentional and harmless here, because nothing in the protocol below depends on the concrete type of the matrix, only on the fact that a frontend can produce one. It is a first hint of the property we will confirm in Section 8.3, namely that an algorithm written once can work with every representation of a [`LatticeModel`](@ref).

## 8.3 Tasks and Data

A frontend produces a momentum-space matrix; what remains is to say *what to compute* from it. That is the task of an assignment, and for the band structure it carries nothing but the reciprocal space to be sampled, quoted from `ToyTBA.jl`:

```julia
struct EigenSystem{R<:ReciprocalSpace}
    reciprocalspace::R
end
function Base.show(io::IO, eigensystem::EigenSystem)
    reciprocalspace = eigensystem.reciprocalspace
    reciprocalspace isa BrillouinZone && return print(io, "EigenSystem(", join(periods(reciprocalspace), "×"), ")")
    return print(io, "EigenSystem")
end
```

Note that the task has no supertype: it is dispatched on as the concrete struct it is, so a marker abstract type would add nothing. The reciprocal space is any subtype of [`ReciprocalSpace`](@ref): a [`BrillouinZone`](@ref) grid, as used below for the density of states, or a [`ReciprocalPath`](@ref), as used for the band-structure plots of [Chapter 7](@ref TutorialAlgorithmInterface). The result of the task is a [`Data`](@ref). For the band structure it stores the eigenvalues as a matrix whose rows are the sampled ``\mathbf{k}``-points and whose columns are the bands, and it keeps the reciprocal space alongside them, so that the result is self-contained for output: the generic recipe of the core package can then plot the assignment directly, as in [Chapter 7](@ref TutorialAlgorithmInterface).

```julia
struct EigenSystemData{R<:ReciprocalSpace} <: Data
    reciprocalspace::R
    values::Matrix{Float64}
end
```

What ties a task to its data is a single method of [`run!`](@ref), dispatched on the pair of an [`Algorithm`](@ref) and an [`Assignment`](@ref):

```julia
function run!(algorithm::Algorithm{<:ToyTBA}, assignment::Assignment{<:EigenSystem}; options...)
    get(options, :showinfo, false) && @info string(assignment)
    bands = Vector{Float64}[]
    for k in assignment.task.reciprocalspace
        push!(bands, eigen(Hermitian(matrix(algorithm, k))).values)
    end
    return EigenSystemData(assignment.task.reciprocalspace, permutedims(reduce(hcat, bands)))
end
```

The method reads what to compute from `assignment.task`; the remaining fields of an [`Assignment`](@ref), namely its parameters, its dependencies and the slot for its data, appear in Section 8.5. The `run!` above is written once, for `Algorithm{<:ToyTBA}`, and yet the frontend it acts on may be either of the two built so far: the method reaches the matrix through `matrix(algorithm, k)` alone, delegated to the frontend as Section 8.2 explained, without ever asking which representation of the model is behind it. The framework, for its part, never asks the developer to declare the result type by hand. Given a task and a frontend, it asks the compiler what `run!` returns:

```@example algorithmdev
datatype(typeof(EigenSystem(path)), typeof(analytic))
```

```@example algorithmdev
datatype(typeof(EigenSystem(path)), typeof(numerical))
```

Both give an `EigenSystemData`, which is the property announced in Section 8.2: one algorithm, one `run!`, one result type, two representations of the model. Two practical consequences follow. First, `run!` must be **type-stable**, otherwise [`datatype`](@ref) cannot infer `Data` and raises an error instead of guessing; moreover, since the inference happens when an assignment is *created*, the error is raised even for a `delay=true` call that computes nothing. Second, this is why [`Data`](@ref) is an abstract type with a prescribed shape rather than a free-form container: it is the contract through which the core package, which knows nothing about any particular algorithm, can still handle its results. In particular, the fields of a `Data` are interpreted in order as the columns of a data file and as the arguments of a plot, so their order matters.

## 8.4 Dependencies

The band structure stands on its own, but the density of states is a sum over eigenvalues: it depends on the eigen-system assignment. Its task, data and `run!` are quoted from `ToyTBA.jl`:

```julia
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
```

The energy window and the broadening are fields of the task itself, with `NaN` as the placeholder for "take it from the data": a `NaN` keeps the fields concretely typed, which a `nothing` default would not. The [`dependencytypes`](@ref) line declares the shape of the dependencies: exactly one, and it must be an assignment of an `EigenSystem` task. The declaration is checked when the assignment is registered on the algorithm, so a malformed call fails at construction time rather than inside `run!`. And the `run!` reads the dependency back from `assignment.dependencies`, so the name used at the call site is never needed inside the implementation. How such a dependency is passed at the call site is part of the assignment machinery of the next section.

## 8.5 The Assignment Machinery

An [`Assignment`](@ref) comes into being when the algorithm is called with a task. The call records the assignment under a name, and since an assignment shares the model interface, it has a directory, an identity and a place to keep its data:

```@example algorithmdev
alg = Algorithm(:Graphene, numerical)
eigensystem = alg(:eigensystem, EigenSystem(BrillouinZone(unitcell, 100)); delay=true)
```

Here the [`BrillouinZone`](@ref) is constructed directly from the lattice, exactly as taught in [Chapter 2](@ref TutorialLattice). This call was made with `delay=true`, so nothing has been computed yet. Calling the algorithm again, now with the assignment, performs the computation:

```@example algorithmdev
alg(eigensystem)
size(eigensystem.data.values)
```

The same call can also be written in the reversed order, `eigensystem(alg)`. The rule that distinguishes the two forms is simple: **the object in parentheses is the authority on parameters**. When the two sides disagree, `alg(assign)` first synchronizes the algorithm, and with it the frontend and the model, to the parameters of the assignment, while `assign(alg)` synchronizes the assignment to those of the algorithm. We use the first form throughout.

Assignments can also depend on one another. The density of states, for instance, takes the eigen-system assignment as a dependency, passed as one more positional argument of the algorithm call:

```@example algorithmdev
dos = alg(:DOS, DensityOfStates(), eigensystem)
```

No parameters were given to this call: an assignment inherits the current parameters of the algorithm (formally, its parameters are `merge(alg.parameters, parameters)`), and per-call values only need to be stated when the task should differ from the algorithm. Even then, the keys of the per-call parameters must be a subset of the parameter keys of the algorithm, so an unknown key is rejected at construction time instead of being silently recorded. This call was made without `delay`, so both the density of states and, if necessary, the eigensystem it depends on have been computed. The result is

```@example algorithmdev
dos.data.energies[1:3]
```

```@example algorithmdev
dos.data.values[1:3]
```

Since no energy range was given, it was taken from the eigenvalues of the dependency. Note also the normalization convention: every eigenvalue contributes a Gaussian of unit area, so the density of states integrates, up to the tails cut off by the energy window, to the total number of sampled single-particle states, here ``2\times100^2``; divide by the number of sampled ``\mathbf{k}`` points if the per-unitcell convention is preferred. Finally, a dependency such as this one is didactic rather than typical: the density of states of the real [TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) package computes the eigenvalues internally and takes no dependency at all.

Underneath these calls sits the caching contract of the interface: an assignment is identified by its task and its parameters, which is what allows the framework to decide by itself whether a result must be computed, can be reused, or has gone stale. What the contract means for the records of a calculation is the subject of the identity and persistence section ([Section 7.5](@ref TutorialPersistence)) of [Chapter 7](@ref TutorialAlgorithmInterface); what it means for the design of a task is taken up in the next section.

## 8.6 Declaring Options

Before the declarations, a distinction that matters. The knobs of a computation come in two kinds. **Definitional** ones change *what* is being computed (the energy window and broadening of a density of states, say) and they live in the task itself, as the fields of `DensityOfStates` showed in Section 8.4. **Execution hints** only change *how* the computation proceeds (tolerances, iteration limits, logging switches) and they travel as keyword arguments. The reason for the split is the caching contract of Section 8.5: an assignment is identified by its task and its parameters, so a definitional knob must be part of the task to be recorded correctly, while hints are free to vary between calls.

The algorithm call accepts keyword arguments, and an algorithm package declares which ones each of its tasks understands. `ToyTBA` declares a single one, quoted from `ToyTBA.jl`:

```julia
@inline options(::Type{<:Assignment{<:EigenSystem}}) = (
    showinfo = "show the information",
)
```

The purpose of the declaration is to make the interface self-describing, and to let the framework police it. A call such as `alg(assign; kwargs...)` is checked against the declared options before anything else happens, so a misspelled keyword is reported together with the list of what was available:

```@example algorithmdev
try
    alg(dos; shwoinfo=true)
catch error
    print(first(error.msg, 55), "...")
end
```

A user who wants the same information before writing the call can simply ask for it. The density of states declares no options of its own, but the options of its dependencies are gathered along the way, which is where the `showinfo` below comes from:

```@example algorithmdev
print(optionsinfo(typeof(dos)))
```

The declaration is not decoration. A `run!` receives the keyword arguments of the call, so `showinfo`, a pure execution hint which changes nothing about what is computed, is honoured by the `run!` of Section 8.3. Two framework rules complete the picture. First, the check can be bypassed by passing `checkoptions=false` to a call, which is exactly what the framework does when it runs the dependencies of a task, so the keyword arguments of a call are in fact forwarded to the `run!` of every assignment in the dependency graph. That is why both `run!` methods of `ToyTBA` end with `options...`: a `run!` must tolerate keywords that are meant for another task. Second, [`update!`](@ref) works on an assignment as well: it replaces the parameters of the task, and the staleness is then detected lazily, by an approximate comparison of the parameter values at the next call. Changing a parameter therefore invalidates the task, so that the next call recomputes both it and, through the algorithm, its dependency:

```@example algorithmdev
update!(dos; t=-1.2)
alg(dos)
dos.data.energies[1:3]
```

Changing a *definitional* knob is a different story. The energy window is part of the task, and the task of an existing assignment is fixed, so a different window means a different assignment. This is the caching contract at work, and it is what keeps every recorded result faithful to the definition that produced it:

```@example algorithmdev
windoweddos = alg(:WindowedDOS, DensityOfStates(emin=-4.0, emax=4.0), eigensystem)
windoweddos.data.energies[1:3]
```

## 8.7 Why the Interface Has This Shape

Now that all the pieces have been written out and exercised, we can look back and see why the interface has the shape it does.

The starting point is a limitation of the language. In Julia only abstract types can be supertypes, and abstract types have no fields. So the state shared by every algorithm, or by every assignment, cannot be hoisted into a common base class; it has to live in **concrete** types. Behaviour, on the other hand, can be inherited: methods defined on an abstract type apply to every subtype. The framework therefore puts the shared *state* into the two concrete types [`Algorithm`](@ref) and [`Assignment`](@ref), and the shared *behaviour* into [`LatticeModel`](@ref), which is why an algorithm and an assignment share the model interface rather than heading hierarchies of their own.

Why, then, do the frontend and the task have no abstract types of their own? Because a marker abstract type enforces nothing: subtyping it guarantees no method, no field, no behaviour. The real contracts of the interface are enforced where they can actually be checked. A frontend must be a [`LatticeModel`](@ref), because an [`Algorithm`](@ref) holds it in a typed field and [`datatype`](@ref) dispatches on it; what it must additionally *do* is fixed by the query methods the algorithm package itself defines and calls. A task declares the shape of its dependencies through [`dependencytypes`](@ref), which is validated when an assignment is registered, and the type of its result is pinned down by [`datatype`](@ref) through the return type of `run!`. Inheritance here is packaging, not taxonomy: a type enters a hierarchy when, and only when, there is shared behaviour to inherit.

| Piece | Kind | State | Why |
|-------|------|-------|-----|
| frontend (role) | a [`LatticeModel`](@ref) subtype per package | not shared | every package's frontend holds something different |
| task (role) | plain struct, no supertype | not shared | every task specifies something different |
| [`Data`](@ref) | abstract | not shared | every result has a different shape |
| [`Algorithm`](@ref) | concrete | shared | every algorithm has the same directory, name, frontend, parameters, `map` and timer |
| [`Assignment`](@ref) | concrete | shared | every assignment has the same directory, name, task, parameters, dependencies and data |

The two concrete types divide the labour between "how to compute" and "what is being computed". An [`Algorithm`](@ref) fixes the method and the parameter values; an [`Assignment`](@ref) fixes a single task, and earns its place by three things an algorithm alone could not provide. It is **named and located**, so its result can be saved, plotted and cited. It keeps its **data in an uninitialised slot**, so that a task can be described without being run, and re-run only when its parameters have changed. And it holds its **dependencies**, which turn a collection of tasks into a computation graph that the framework walks in the right order.

That leaves the question of why [`Data`](@ref) is designed the way it is, and the answer is subtler than it first appears.

The decisive reason is that an [`Assignment`](@ref) is constructed **before** anything is computed, and may never be computed at all if `delay=true`. Its result slot must nevertheless have a concrete type, and Julia allows a field to be left undefined only if its type is known at construction time. The type of the result therefore cannot be decided after the fact; it has to be derived statically from the pair of the frontend and the task, which is exactly what [`datatype`](@ref) does by asking the compiler for the return type of [`run!`](@ref). This is why `Data` must be an abstract type with a well-defined interface, why the type of the result is a type parameter of the [`Assignment`](@ref) rather than a property of the task, and why `run!` is required to be type-stable.

The second reason is genericity of the operations that come after the computation. The `Tuple` of a `Data` returns its fields in order, and both the plain-text output and the plotting recipes are defined in terms of that tuple. The core package can therefore save and draw the results of an algorithm it has never heard of, in the same way that it can store and stamp a model it has never heard of. Note the asymmetry with the task: a `Data` *does* have shared behaviour worth packaging (the tuple view and equality), so it gets an abstract type; a task has none, so it gets none.

The third reason is that a result is not a description. A `Data` has no parameters, no configuration, and no identity, so the model interface would be largely dead weight on it. It keeps only what is needed: equality for comparison, and the tuple view for output.

## Summary

The contract of an algorithm package, exercised in full on `ToyTBA`:

- A frontend, a [`LatticeModel`](@ref) subtype defined by the algorithm package, encapsulates what one algorithm needs from a model, and because it is written against [`LatticeModel`](@ref), a single algorithm works with every representation of a model; [`@delegate`](@ref) mirrors its query interface on the algorithm.
- A task and a [`Data`](@ref) say what to compute and what comes out, and the two are wired together automatically, [`datatype`](@ref) inferring the latter from the return type of `run!`, which must therefore be type-stable, while [`dependencytypes`](@ref) declares the shape of the dependencies and is checked when an assignment is registered.
- [`options`](@ref) and [`optionsinfo`](@ref) make the keyword arguments of a call self-documenting and checkable; options are execution hints and never take part in the identity of an assignment, while definitional knobs belong to the task itself.

The complete source of the toy package dissected in this chapter is the file `ToyTBA.jl` included in this documentation, and the docstrings of the types and functions introduced in this chapter are collected in the [manual](@ref "Frameworks").
