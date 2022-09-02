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
export dtype, kind, reset!, update, update!

# Essentials.QuantumOperators
export ID, Operator, OperatorPack, OperatorProd, Operators, OperatorSum, OperatorUnit, LaTeX, QuantumOperator
export AbstractSubstitution, AbstractUnitSubstitution, Identity, MatrixRepresentation, Numericalization, Permutation, RankFilter, Transformation, UnitSubstitution
export idtype, ishermitian, latexname, latexformat, matrix, matrix!, optype, script, sequence, subscript, superscript

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, Momentum, Momentum₁, Momentum₂, Momentum₃, ParticleNumber, SpinfulParticle, SpinZ
export particlenumbers, periods, spinfulparticles, spinzs, @abeliannumber

# Essentials.Spatials
export azimuth, azimuthd, distance, isintratriangle, isonline, isparallel, issubordinate, interlinks, minimumlengths, polar, polard, reciprocals, rotate, translate, tile, volume
export AbstractLattice, Bond, BrillouinZone, Lattice, Neighbors, Point, ReciprocalSpace, ReciprocalZone, ReciprocalPath, Segment, Translations, bonds!, bonds, icoordinate, isintracell, nneighbor, rcoordinate, @translations_str
export hexagon120°map, hexagon60°map, linemap, rectanglemap, @hexagon_str, @line_str, @rectangle_str

# Essentials.DegreesOfFreedom
export Boundary, CompositeIID, CompositeIndex, CompositeInternal, Coupling, Couplings, Hilbert, IID, IIDSpace, Index, Internal, Metric, OperatorUnitToTuple, SimpleIID, SimpleInternal, Subscript, Subscripts, Table, Term
export plain, iidtype, indextype, ismodulatable, statistics, @couplings, @subscript_str

# Essentials.QuantumSystems
## Canonical fermionic/bosonic systems
export annihilation, creation, latexofbosons, latexoffermions, majorana
export Coulomb, FID, Fock, FockCoupling, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip
export isnormalordered, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str

## SU(2) spin systems
export latexofspins
export SID, Spin, SpinCoupling, SpinTerm
export totalspin, @dm_str, @gamma_str, @heisenberg_str, @ising_str, @sc_str, @sˣ_str, @sʸ_str, @sᶻ_str

## Phononic systems
export latexofphonons
export PID, Phonon, PhononCoupling, PhononKinetic, PhononPotential, PhononTerm
export @kinetic_str, @potential_str

# Essentials.Frameworks
export Action, Algorithm, AnalyticalExpression, Assignment, CompositeGenerator, Entry, Frontend, Image, OperatorGenerator, Parameters, RepresentationGenerator, prepare!, run!, rundependences!, save

end
