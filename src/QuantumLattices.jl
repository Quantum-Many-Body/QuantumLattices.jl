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
export QuantumOperator, OperatorUnit, ID, OperatorPack, OperatorProd, OperatorSum, Operator, Operators, LaTeX
export Transformation, Identity, Numericalization, MatrixRepresentation, Permutation, AbstractSubstitution, AbstractUnitSubstitution, UnitSubstitution, RankFilter
export ishermitian, idtype, optype, sequence, latexname, latexformat, matrix, matrix!

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, @abeliannumber, periods, Momentum, Momentum₁, Momentum₂, Momentum₃

# Essentials.Spatials
export distance, azimuthd, azimuth, polard, polar, volume, isparallel, isonline, isintratriangle, issubordinate, reciprocals, translate, rotate, tile
export Translations, AbstractPID, PID, CPID, AbstractBond, Point, Bond, Lattice, SuperLattice, Cylinder, Bonds
export Segment, ReciprocalSpace, BrillouinZone, ReciprocalZone, ReciprocalPath
export pidtype, rcoord, icoord, isintracell, bonds!, bonds, @line_str, @rectangle_str, @hexagon_str, @translations_str
export allbonds, zerothbonds, insidebonds, acrossbonds, intrabonds, interbonds

# Essentials.DegreesOfFreedom
export IID, SimpleIID, CompositeIID, Internal, SimpleInternal, CompositeInternal, AbstractOID, Index, CompositeOID, OID, Hilbert
export IIDSpace, Subscript, Subscripts, AbstractCoupling, Coupling, Couplings, Metric, OIDToTuple, Table, Term
export statistics, iidtype, indextype, ismodulatable, abbr, @subscript_str, @couplings

# Essentials.Frameworks
export Parameters, Boundary, Engine, AbstractGenerator, Formulation, Entry, CompositeGenerator, Generator, Image, Action, Assignment, Algorithm
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

end
