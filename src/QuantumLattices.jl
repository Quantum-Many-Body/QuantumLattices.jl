module QuantumLattices

include("Interfaces.jl")
include("Prerequisites/Prerequisites.jl")
include("Essentials/Essentials.jl")

using .Interfaces
using .Prerequisites
using .Essentials

# Interfaces
export dimension, id, rank, value
export ⊕, ⊗, ⋅, add!, div!, mul!, sub!
export decompose, decompose!, expand, expand!, permute

# Essentials
export kind, reset!, update, update!

# Essentials.QuantumOperators
export ID, Operator, OperatorPack, OperatorProd, Operators, OperatorSum, OperatorUnit, LaTeX, QuantumOperator
export AbstractSubstitution, AbstractUnitSubstitution, Identity, MatrixRepresentation, Numericalization, Permutation, RankFilter, Transformation, UnitSubstitution
export ishermitian, latexname, latexformat, matrix, matrix!, script, sequence

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, Momentum, Momentum₁, Momentum₂, Momentum₃, ParticleNumber, SpinfulParticle, SpinZ
export particlenumbers, periods, spinfulparticles, spinzs, @abeliannumber

# Essentials.Spatials
export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, Segment, Translations, bonds!, bonds, icoordinate, isintracell, rcoordinate, @translations_str
export @hexagon_str, @line_str, @rectangle_str

# Essentials.DegreesOfFreedom
export Boundary, CompositeIID, CompositeIndex, CompositeInternal, Component, Constraint, Coupling, Hilbert, IID, IIDSpace, Index, Internal, MatrixCoupling, Metric, OperatorUnitToTuple, SimpleIID, SimpleInternal, Table, Term
export plain, statistics, @indexes

# Essentials.QuantumSystems
## Canonical complex fermionic/bosonic systems
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, annihilation, creation, latexofbosons, latexoffermions, latexofparticles
export Coulomb, FID, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isnormalordered

## SU(2) spin systems
export latexofspins, SID, Spin, SpinTerm, totalspin, @dm_str, @gamma_str, @heisenberg_str, @ising_str

## Phononic systems
export latexofphonons, Elastic, PID, Phonon, Kinetic, Hooke, PhononTerm

# Essentials.Frameworks
export Action, Algorithm, AnalyticalExpression, Assignment, CompositeGenerator, Entry, Frontend, Image, OperatorGenerator, Parameters, RepresentationGenerator, prepare!, run!, rundependences!, save

end
