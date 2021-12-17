module QuantumLattices

include("Interfaces.jl")
include("Prerequisites/Prerequisites.jl")
include("Essentials/Essentials.jl")

using .Interfaces
using .Prerequisites
using .Essentials

# Interfaces
export id, value, rank, dimension
export ⊕, ⊗, ⋅, add!, sub!, mul!, div!
export expand, expand!, decompose, decompose!, permute

# Essentials
export dtype, kind, update, update!, reset!

# Essentials.QuantumOperators
export QuantumOperator, OperatorUnit, ID, OperatorProd, OperatorSum, Operator, Operators, ishermitian, idtype, optype, sequence, latexformat
export Transformation, Identity, Numericalization, MatrixRepresentation, Permutation, UnitSubstitution, RankFilter, matrix, matrix!

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, @abeliannumber, periods, Momentum, Momentum₁, Momentum₂, Momentum₃

# Essentials.Spatials
export distance, azimuthd, azimuth, polard, polar, volume, isparallel, isonline, isintratriangle, issubordinate, reciprocals, translate, rotate
export Translations, PID, CPID, Point, Bond, Lattice, SuperLattice, Cylinder, Bonds, Segment, BrillouinZone, ReciprocalZone, ReciprocalPath
export pidtype, rcoord, icoord, isintracell, bonds!, bonds, @line_str, @rectangle_str, @hexagon_str, @translations_str
export allbonds, zerothbonds, insidebonds, acrossbonds, intrabonds, interbonds

# Essentials.DegreesOfFreedom
export CompositeIID, CompositeInternal, Index, OID, IIDSpace, Hilbert
export Subscript, Subscripts, SubscriptsID, Coupling, Couplings, OIDToTuple, Table, Term, LaTeX
export statistics, iidtype, ismodulatable, abbr, @subscript_str, @couplings

# Essentials.Frameworks
export Parameters, Boundary, Engine, Formulation, Entry, Generator, Image, Action, Assignment, Algorithm
export prepare!, run!, rundependences!, plain

# Essentials.QuantumSystems
## Canonical fermionic/bosonic systems
export majorana, annihilation, creation, flatex, blatex
export FID, Fock, FockCoupling, Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb, FockTerm
export isnormalordered, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str

## SU(2) spin systems
export slatex
export SID, Spin, SpinCoupling, SpinTerm, totalspin
export @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

## Phononic systems
export nlatex
export NID, Phonon, PhononCoupling, PhononKinetic, PhononPotential, PhononTerm

## Magnon-Phonon coupled systems
export DMPhonon

end
