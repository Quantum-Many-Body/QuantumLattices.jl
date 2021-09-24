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
export dtype, kind, update!, reset!

# Essentials.QuantumOperators
export ID, QuantumOperator, OperatorProd, OperatorSum, Transformation, Numericalization, MatrixRepresentation, idtype, sequence, matrix

# Essentials.QuantumNumbers
export AbelianNumber, AbelianNumbers, @abeliannumber, periods, Momentum, Momentum₁, Momentum₂, Momentum₃

# Essentials.Spatials
export distance, azimuthd, azimuth, polard, polar, volume, isparallel, isonline, isintratriangle, issubordinate, reciprocals, translate, rotate
export PID, CPID, Point, Bond, Lattice, SuperLattice, Cylinder, Bonds, Segment, BrillouinZone, ReciprocalZone, ReciprocalPath
export pidtype, rcoord, icoord, isintracell, bonds!, bonds, @line_str, @rectangle_str, @hexagon_str
export allbonds, zerothbonds, insidebonds, acrossbonds, intrabonds, interbonds

# Essentials.DegreesOfFreedom
export CompositeIID, CompositeInternal, Index, OID, Operator, Operators, IIDSpace, Hilbert
export Subscript, Subscripts, SubscriptsID, Coupling, Couplings, OIDToTuple, Table, Term, LaTeX, Boundary
export statistics, iidtype, ishermitian, ismodulatable, abbr, otype, latexformat, twist, plain, @subscript_str, @couplings

# Essentials.Frameworks
export Parameters, Generator, SimplifiedGenerator, Action, Engine, Assignment, Algorithm
export prepare!, register!, run!, dependences, rundependences!

# Essentials.QuantumSystems
## Canonical fermionic/bosonic systems
export majorana, annihilation, creation, flatex, blatex
export FID, Fock, FockCoupling, Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb
export isnormalordered, @σ⁰_str, @σˣ_str, @σʸ_str, @σᶻ_str, @σ⁺_str, @σ⁻_str, @fc_str

## SU(2) spin systems
export slatex
export SID, Spin, SpinCoupling, SpinTerm, totalspin
export @heisenberg_str, @ising_str, @gamma_str, @dm_str, @sˣ_str, @sʸ_str, @sᶻ_str, @sc_str

## Phononic systems
export nlatex
export NID, Phonon, PhononCoupling, PhononKinetic, PhononPotential

## Magnon-Phonon coupled systems
export DMPhonon

end
