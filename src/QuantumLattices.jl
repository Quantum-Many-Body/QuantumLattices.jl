module QuantumLattices

# interfaces
import LinearAlgebra: ishermitian, mul!, rank
include("Interfaces.jl")
export OneAtLeast, OneOrMore, ZeroAtLeast, ZeroOrMore
export ⊕, ⊗, add!, decompose, decompose!, dimension, div!, expand, expand!, id, ishermitian, kind, mul!, permute, rank, reset!, shape, str, sub!, update, update!, value

# Toolkit
include("Toolkit.jl")
using .Toolkit

# QuantumOperators
include("QuantumOperators.jl")
using .QuantumOperators
export LaTeX, Operator, OperatorIndex, OperatorPack, OperatorProd, Operators, OperatorSet, OperatorSum, QuantumOperator
export LinearFunction, LinearTransformation, Matrixization, Permutation, TabledUnitSubstitution, UnitSubstitution
export idtype, isequivalenttoscalar, ishermitian, latexname, latexformat, matrix, scalartype, script, sequence

# Spatials
include("Spatials.jl")
using .Spatials
export azimuth, azimuthd, direction, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, nneighbor, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, FractionalReciprocalSpace, Lattice, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalSpace, ReciprocalZone
export bonds, bonds!, dlmsave, fractionals, icoordinate, iscontinuous, isdiscrete, isintracell, label, period, periods, rcoordinate, selectpath, shrink, ticks, @hexagon_str, @line_str, @rectangle_str

# DegreesOfFreedom
include("DegreesOfFreedom.jl")
using .DegreesOfFreedom
export CompositeInternal, CompositeIndex, CoordinatedIndex, Hilbert, Index, Internal, InternalIndex, InternalProd, InternalSum, SimpleInternal
export Boundary, Coupling, MatrixCoupling, MatrixCouplingComponent, MatrixCouplingProd, MatrixCouplingSum, Metric, OperatorIndexToTuple, Ordinal, Pattern, Table, Term, TermAmplitude, TermCoupling
export ˢᵗ, ⁿᵈ, ʳᵈ, ᵗʰ, plain, coordinatedindextype, diagonalfields, indextype, internalindextype, isdefinite, isdiagonal, partition, patternrule, showablefields, statistics, @pattern

# QuantumSystems
include("QuantumSystems.jl")
using .QuantumSystems
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, σ¹¹, σ¹², σ²¹, σ²², annihilation, creation, latexofbosons, latexoffermions, latexofparticles, Lˣ, Lʸ, Lᶻ
export 𝕒, 𝕒⁺, 𝕒𝕒, 𝕒𝕒⁺, 𝕒⁺𝕒, 𝕒⁺𝕒⁺, 𝕔, 𝕔⁺, 𝕔𝕔, 𝕔𝕔⁺, 𝕔⁺𝕔, 𝕔⁺𝕔⁺, 𝕕, 𝕕⁺, 𝕕𝕕, 𝕕𝕕⁺, 𝕕⁺𝕕, 𝕕⁺𝕕⁺, Fock, FockIndex, isannihilation, iscreation, isnormalordered
export Coulomb, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip
export Γˣ, Γʸ, Γᶻ, Γ′ˣ, Γ′ʸ, Γ′ᶻ, DMˣ, DMʸ, DMᶻ, Isingˣ, Isingʸ, Isingᶻ, latexofspins
export 𝕊, 𝕊ᵀ𝕊, SpinIndex, Spin, totalspin
export Γ, Γ′, DM, Heisenberg, Ising, Kitaev, SingleIonAnisotropy, SpinTerm, Zeeman
export latexofphonons
export 𝕦, 𝕦ᵀ𝕦, 𝕡, Phonon, PhononIndex
export Elastic, Kinetic, Hooke, PhononTerm

# Frameworks
include("Frameworks.jl")
using .Frameworks
export Action, Algorithm, Assignment, CategorizedGenerator, Data, Eager, ExpansionStyle, Formula, Frontend, Generator, Lazy, OperatorGenerator, Parameters
export eager, lazy, checkoptions, datatype, hasoption, options, optionsinfo, qldload, qldsave, run!

end
