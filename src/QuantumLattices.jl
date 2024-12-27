module QuantumLattices

import LinearAlgebra: ishermitian, rank, mul!

# interfaces
export ⊞, ⊠, ⊕, ⊗, add!, decompose, decompose!, dimension, div!, expand, expand!, id, ishermitian, kind, mul!, permute, rank, reset!, sub!, update, update!, value

function ⊞ end
function ⊠ end
function ⊕ end
function ⊗ end
function add! end
function decompose end
function decompose! end
function dimension end
function div! end
function expand end
function expand! end
function id end
function kind end
function permute end
function reset! end
function sub! end
function update end
function update! end
function value end

# Toolkit
include("Toolkit.jl")
using .Toolkit

# QuantumOperators
include("QuantumOperators.jl")
using .QuantumOperators
export ID, LaTeX, Operator, OperatorIndex, OperatorPack, OperatorProd, Operators, OperatorSet, OperatorSum, QuantumOperator
export LinearFunction, LinearTransformation, Matrixization, Permutation, RankFilter, TabledUnitSubstitution, Transformation, UnitSubstitution
export equivalenttoscalar, idtype, ishermitian, latexname, latexformat, matrix, scalartype, script, sequence

# QuantumNumbers
include("QuantumNumbers.jl")
using .QuantumNumbers
export Abelian, AbelianQuantumNumber, AbelianGradedSpace, AbelianGradedSpaceProd, AbelianGradedSpaceSum, CompositeAbelianQuantumNumber, Graded, Momenta, RepresentationSpace, SimpleAbelianQuantumNumber
export 𝕂, 𝕂¹, 𝕂², 𝕂³, ℕ, 𝕊ᶻ, 𝕌₁, ℤ, ℤ₁, ℤ₂, ℤ₃, ℤ₄, findindex, period, periods, regularize, regularize!

# Spatials
include("Spatials.jl")
using .Spatials
export azimuth, azimuthd, direction, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, nneighbor, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalCurve, ReciprocalSpace, ReciprocalZone, ReciprocalPath, bonds, bonds!, icoordinate, isintracell, rcoordinate, save, selectpath, shrink
export @hexagon_str, @line_str, @rectangle_str

# DegreesOfFreedom
include("DegreesOfFreedom.jl")
using .DegreesOfFreedom
export AllEqual, Boundary, CompositeInternal, ConstrainedInternal, Internal, InternalIndex, InternalIndexProd, InternalPattern, InternalProd, InternalSum, SimpleInternal, SimpleInternalIndex
export Component, CompositeIndex, CoordinatedIndex, Coupling, Hilbert, Index, MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum, Metric, OperatorIndexToTuple, Ordinal, Pattern, Table, Term, TermAmplitude, TermCoupling, TermFunction
export ˢᵗ, ⁿᵈ, ʳᵈ, ᵗʰ, plain, allequalfields, coordinatedindextype, indextype, internalindextype, isdefinite, partition, patternrule, sitestructure, statistics, @pattern

# QuantumSystems
include("QuantumSystems.jl")
using .QuantumSystems
export annihilation, creation, latexofbosons, latexoffermions, latexofparticles, 𝕒, 𝕒⁺𝕒, 𝕓, 𝕓⁺𝕓, 𝕔, 𝕔⁺𝕔, 𝕕, 𝕕⁺𝕕, 𝕗, 𝕗⁺𝕗, isannihilation, iscreation, isnormalordered, @σ_str, @L_str
export Coulomb, Fock, FockIndex, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip
export latexofspins, 𝕊, 𝕊ᵀ𝕊, SpinIndex, Spin, totalspin, @Γ_str, @Γ′_str, @DM_str, @Heisenberg_str, @Ising_str
export DM, Heisenberg, Ising, Kitaev, SingleIonAnisotropy, SpinTerm, Zeeman, Γ, Γ′
export latexofphonons, Elastic, Phonon, PhononIndex, Kinetic, Hooke, PhononTerm, 𝕦, 𝕦ᵀ𝕦, 𝕡

# Frameworks
include("Frameworks.jl")
using .Frameworks
export Action, Algorithm, Assignment, CategorizedGenerator, Eager, ExpansionStyle, Formula, Frontend, Generator, Lazy, OneOrMore, OperatorGenerator, Parameters
export checkoptions, eager, lazy, initialize, options, prepare!, run!

end
